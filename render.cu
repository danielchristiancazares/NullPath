#define BLACKHOLE_NO_MAIN
#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>

#include "bh_common.h"
#include "bh_types.h"
#include "bh_schwarzschild.h"
#include "bh_kerr.h"

struct RenderParams {
    int width = 1280;
    int height = 720;
    float fov_y_deg = 50.0f;
    float theta_deg = 90.0f; // camera polar angle
    float mass = 10.0f;     // solar masses in geometric units (M)
    float r_cam = 50.0f;    // camera radius in units of Rs
    const char* out_path = "render.ppm";
    int samples = 1;        // samples per pixel
    int tile_h = 128;       // rows per tile for async copy/compute
    int frames = 1;         // number of frames
    float turns = 1.0f;     // camera turns around BH over all frames
    float phi0_deg = 0.0f;  // starting azimuth
    const char* out_pat = "frame_%05d.ppm"; // used when frames > 1
    float spin = 0.0f;      // a/M in [0,1)
    int spp_total = 0;      // total spp across passes (0 => just 'samples')
    int checkpoint_every = 0; // checkpoint every K passes (0 => final only)
    int check_invariants = 0; // 1 = track Hamiltonian/E/Lz drift and report
    const char* mode = "balanced"; // fast|balanced|accurate
};

#ifndef BH_FP64_GEODESIC
#define BH_FP64_GEODESIC 1
#endif

__device__ __host__ inline float3 make_float3_clamped(float x, float y, float z) {
    return make_float3(fminf(fmaxf(x, 0.f), 1.f), fminf(fmaxf(y, 0.f), 1.f), fminf(fmaxf(z, 0.f), 1.f));
}

// (Kerr metric, camera init, and Doppler are provided by headers)

__device__ inline float3 shade_disk(float r, float Rs, float g) {
    // Emissivity j(r) ~ r^{-2}; apply relativistic beaming: I_obs = g^3 I_em
    float j = powf(fmaxf(r, 1.0f), -2.0f);
    float I = j * g * g * g;
    float3 base = make_float3(1.0f, 0.75f, 0.45f);
    return make_float3_clamped(base.x * I, base.y * I, base.z * I);
}

__device__ inline float3 shade_background(float u, float v) {
    // Simple sky gradient with a bright band to mimic stars
    float t = fminf(fmaxf(v, 0.f), 1.f);
    float3 top = make_float3(0.02f, 0.03f, 0.05f);
    float3 bottom = make_float3(0.08f, 0.09f, 0.12f);
    float3 c = make_float3(top.x*(1-t) + bottom.x*t,
                           top.y*(1-t) + bottom.y*t,
                           top.z*(1-t) + bottom.z*t);
    // Add faint band
    float band = expf(-50.0f * (v-0.5f)*(v-0.5f));
    c.x += 0.4f * band; c.y += 0.35f * band; c.z += 0.3f * band;
    return make_float3_clamped(c.x, c.y, c.z);
}

__device__ inline unsigned int wanghash(unsigned int a) {
    a = (a ^ 61u) ^ (a >> 16);
    a *= 9u;
    a = a ^ (a >> 4);
    a *= 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}

__device__ inline float rnd01(unsigned int seed) {
    return (wanghash(seed) & 0x00ffffff) * (1.0f / 16777216.0f);
}

// Atomic max for float via CAS
__device__ inline void atomicMaxFloat(float* addr, float val) {
    int* address_as_i = (int*)addr;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        float old_f = __int_as_float(assumed);
        float max_f = fmaxf(old_f, val);
        old = atomicCAS(address_as_i, assumed, __float_as_int(max_f));
    } while (assumed != old);
}

// (Hamiltonians provided by headers)

// (Kerr derivatives/integrator provided by headers)

// Kerr-aware render kernel
__global__ void renderKernel(int width, int height, int y0, int tile_h, int spp,
                             float mass, float spin_a_over_M, float r_cam_Rs, float theta_cam, float phi_cam,
                             float fov_y_deg, int check_inv, int integrator_mode, int camera_mode,
                             float* inv_max3, uchar3* out_tile) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tileN = width * tile_h;
    if (idx >= tileN) return;
    int x = idx % width;
    int y = y0 + (idx / width);
    if (y >= height) return;

    SchwarzschildMetric metric(mass);
    float Rs = metric.Rs; // still equals 2M
    float M = mass;
    float a = spin_a_over_M * M;
    float r_cam = r_cam_Rs * Rs;

    float aspect = (float)width / (float)height;
    float fov_y = fov_y_deg * BH_PI / 180.0f;
    float tan_half = tanf(0.5f * fov_y);

    float3 accum = make_float3(0,0,0);
    for (int s = 0; s < spp; ++s) {
        // Subpixel jitter
        float jx = rnd01((unsigned)(x*73856093u ^ y*19349663u ^ s*83492791u)) - 0.5f;
        float jy = rnd01((unsigned)(x*83492791u ^ y*2654435761u ^ s*19349663u)) - 0.5f;

        float ndc_x = ((x + 0.5f + jx) / (float)width) * 2.0f - 1.0f;
        float ndc_y = ((y + 0.5f + jy) / (float)height) * 2.0f - 1.0f;
        float cam_x = ndc_x * aspect * tan_half;
        float cam_y = -ndc_y * tan_half;
        float cam_z = -1.0f;
        float norm = rsqrtf(cam_x*cam_x + cam_y*cam_y + cam_z*cam_z);
        cam_x *= norm; cam_y *= norm; cam_z *= norm;

        // Camera basis
        float n_r = -cam_z;
        float n_theta = cam_y;
        float n_phi = cam_x;

        #if BH_FP64_GEODESIC
        GeodesicStateD state;
        bool static_ok = kerr_static_timelike((double)M, (double)a, (double)r_cam, (double)theta_cam);
        if (camera_mode == 1 || !static_ok) {
            state = init_ray_from_camera_zamo_d((double)M, (double)a, (double)r_cam, (double)theta_cam, (double)phi_cam,
                                                (double)n_r, (double)n_theta, (double)n_phi);
        } else {
            state = init_ray_from_camera_general_d((double)M, (double)a, (double)r_cam, (double)theta_cam, (double)phi_cam,
                                                   (double)n_r, (double)n_theta, (double)n_phi);
        }
        double pt0 = state.pt;
        double pphi0 = state.pphi;
        double min_r = state.r;
        bool captured = false;
        bool disk_hit = false;
        double hit_r = 0.0;
        double prev_theta = state.theta;
        GeodesicStateD prev_state = state;
        const double rin = 3.0 * (double)Rs;
        const double rout = 20.0 * (double)Rs;
        float g_boost = 1.0f;
        for (int step = 0; step < MAX_INTEGRATION_STEPS; ++step) {
            double scale = fmin(fmax(1.0, state.r / (0.5 * (double)Rs)), 100.0);
            bool ok = false;
            if (integrator_mode == 1) {
                double h = (double)INTEGRATION_STEP_SIZE * scale;
                h = fmin(fmax(h, BH_RK_HMIN), BH_RK_HMAX);
                // Perform one accepted adaptive step (may shrink h internally)
                double h_try = h;
                double fac = rk4_adaptive_step_d(state, (double)M, (double)a, h_try, BH_RK_REL_TOL, BH_RK_ABS_TOL);
                ok = (fac >= 0.0); // accepted
            } else {
                ok = integrateGeodesicKerrD(state, (double)M, (double)a, (double)INTEGRATION_STEP_SIZE * scale);
            }
            if (!ok) { captured = true; break; }
            min_r = fmin(min_r, state.r);
            double d_prev = prev_state.theta - (double)BH_PI_D*0.5;
            double d_curr = state.theta - (double)BH_PI_D*0.5;
            if (!disk_hit) {
                double eps = 1e-8;
                bool sign_change = (d_prev > 0 && d_curr < 0) || (d_prev < 0 && d_curr > 0);
                bool near_zero_transition = (fabs(d_prev) < eps && fabs(d_curr) >= eps) || (fabs(d_curr) < eps && fabs(d_prev) >= eps);
                if (sign_change || near_zero_transition) {
                    double denom = fabs(d_prev) + fabs(d_curr) + 1e-30;
                    double alpha = fabs(d_prev) / denom;
                    double r_cross = prev_state.r + alpha * (state.r - prev_state.r);
                    if (r_cross >= rin && r_cross <= rout) {
                        GeodesicStateD cross = state;
                        cross.r = r_cross;
                        cross.theta = (double)BH_PI_D*0.5;
                        cross.pt = prev_state.pt + alpha * (state.pt - prev_state.pt);
                        cross.pphi = prev_state.pphi + alpha * (state.pphi - prev_state.pphi);
                        disk_hit = true;
                        hit_r = r_cross;
                        g_boost = (float)doppler_g_general_d(cross, (double)M, (double)a);
                        if (check_inv) {
                            // Recompute H in double
                            double gtt, gtphi, gphiphi, grr, gthth;
                            kerr_blocks_contra_d((double)M, (double)a, cross.r, cross.theta, gtt, gtphi, gphiphi, grr, gthth);
                            double Hc = 0.5*(gtt*cross.pt*cross.pt + 2.0*gtphi*cross.pt*cross.pphi + gphiphi*cross.pphi*cross.pphi
                                             + grr*cross.pr*cross.pr + gthth*cross.ptheta*cross.ptheta);
                            float Habs = (float)fabs(Hc);
                            float Erel = (float)(fabs(cross.pt - pt0) / fmax(1e-12, fabs(pt0)));
                            float Lrel = (float)(fabs(cross.pphi - pphi0) / fmax(1e-12, fabs(pphi0)));
                            atomicMaxFloat(&inv_max3[0], Habs);
                            atomicMaxFloat(&inv_max3[1], Erel);
                            atomicMaxFloat(&inv_max3[2], Lrel);
                        }
                    }
                }
            }
            prev_theta = state.theta;
            prev_state = state;
            if (state.r > 999.0 && state.pr > 0) { break; }
        }
        #else
        GeodesicState state = init_ray_from_camera_general(M, a, r_cam, theta_cam, phi_cam, n_r, n_theta, n_phi);
        float pt0 = state.pt;
        float pphi0 = state.pphi;
        float min_r = state.r;
        bool captured = false;
        bool disk_hit = false;
        float hit_r = 0.0f;
        float prev_theta = state.theta;
        GeodesicState prev_state = state;
        const float rin = 3.0f * Rs;
        const float rout = 20.0f * Rs;
        float g_boost = 1.0f;
        for (int step = 0; step < MAX_INTEGRATION_STEPS; ++step) {
            bool ok = false;
            float scale = fminf(100.0f, fmaxf(1.0f, state.r / (0.5f * Rs)));
            if (a == 0.0f) ok = integrateGeodesic(state, metric, INTEGRATION_STEP_SIZE * scale);
            else ok = integrateGeodesicKerr(state, M, a, INTEGRATION_STEP_SIZE * scale);
            if (!ok) { captured = true; break; }
            min_r = fminf(min_r, state.r);
            float d_prev = prev_state.theta - BH_PI*0.5f;
            float d_curr = state.theta - BH_PI*0.5f;
            if (!disk_hit) {
                float eps = 1e-4f;
                bool sign_change = (d_prev > 0 && d_curr < 0) || (d_prev < 0 && d_curr > 0);
                bool near_zero_transition = (fabsf(d_prev) < eps && fabsf(d_curr) >= eps) || (fabsf(d_curr) < eps && fabsf(d_prev) >= eps);
                if (sign_change || near_zero_transition) {
                    float denom = fabsf(d_prev) + fabsf(d_curr) + 1e-20f;
                    float alpha = fabsf(d_prev) / denom;
                    float r_cross = prev_state.r + alpha * (state.r - prev_state.r);
                    if (r_cross >= rin && r_cross <= rout) {
                        GeodesicState cross = state;
                        cross.r = r_cross; cross.theta = PI_F*0.5f;
                        cross.pt = prev_state.pt + alpha * (state.pt - prev_state.pt);
                        cross.pphi = prev_state.pphi + alpha * (state.pphi - prev_state.pphi);
                        disk_hit = true; hit_r = r_cross;
                        g_boost = doppler_g_general(cross, M, a);
                        if (check_inv) {
                            float Hc = (a == 0.0f) ? hamiltonian_schwarzschild(cross, Rs) : hamiltonian_kerr(cross, M, a);
                            float Habs = fabsf(Hc);
                            float Erel = fabsf(cross.pt - pt0) / fmaxf(1e-9f, fabsf(pt0));
                            float Lrel = fabsf(cross.pphi - pphi0) / fmaxf(1e-9f, fabsf(pphi0));
                            atomicMaxFloat(&inv_max3[0], Habs);
                            atomicMaxFloat(&inv_max3[1], Erel);
                            atomicMaxFloat(&inv_max3[2], Lrel);
                        }
                    }
                }
            }
            prev_theta = state.theta;
            prev_state = state;
            if (state.r > 999.0f && state.pr > 0) { break; }
        }
        #endif

        float3 color;
        if (captured) {
            color = make_float3(0,0,0);
        } else if (disk_hit) {
            color = shade_disk((float)hit_r, Rs, g_boost);
        } else {
            if (check_inv) {
                #if BH_FP64_GEODESIC
                double gtt, gtphi, gphiphi, grr, gthth;
                kerr_blocks_contra_d((double)M, (double)a, state.r, state.theta, gtt, gtphi, gphiphi, grr, gthth);
                double Hf = 0.5*(gtt*state.pt*state.pt + 2.0*gtphi*state.pt*state.pphi + gphiphi*state.pphi*state.pphi
                                  + grr*state.pr*state.pr + gthth*state.ptheta*state.ptheta);
                float Habs = (float)fabs(Hf);
                float Erel = (float)(fabs(state.pt - pt0) / fmax(1e-12, fabs(pt0)));
                float Lrel = (float)(fabs(state.pphi - pphi0) / fmax(1e-12, fabs(pphi0)));
                #else
                float Hf = (a == 0.0f) ? hamiltonian_schwarzschild(state, Rs) : hamiltonian_kerr(state, M, a);
                float Habs = fabsf(Hf);
                float Erel = fabsf(state.pt - pt0) / fmaxf(1e-9f, fabsf(pt0));
                float Lrel = fabsf(state.pphi - pphi0) / fmaxf(1e-9f, fabsf(pphi0));
                #endif
                atomicMaxFloat(&inv_max3[0], Habs);
                atomicMaxFloat(&inv_max3[1], Erel);
                atomicMaxFloat(&inv_max3[2], Lrel);
            }
            float u = (cam_x * 0.5f + 0.5f);
            float v = (cam_y * 0.5f + 0.5f);
            color = shade_background(u, v);
        }

        accum.x += color.x; accum.y += color.y; accum.z += color.z;
    }
    float inv = 1.0f / (float)spp;
    float3 color = make_float3(accum.x * inv, accum.y * inv, accum.z * inv);
    out_tile[idx] = make_uchar3((unsigned char)(255.0f * color.x),
                                (unsigned char)(255.0f * color.y),
                                (unsigned char)(255.0f * color.z));
}

static bool write_ppm(const char* path, int w, int h, const std::vector<unsigned char>& rgb) {
    FILE* f = fopen(path, "wb");
    if (!f) return false;
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    size_t wrote = fwrite(rgb.data(), 1, rgb.size(), f);
    fclose(f);
    return wrote == rgb.size();
}

static void parse_args(int argc, char** argv, RenderParams& p) {
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        auto nextf = [&](float& ref){ if (i+1 < argc) ref = strtof(argv[++i], nullptr); };
        auto nexti = [&](int& ref){ if (i+1 < argc) ref = atoi(argv[++i]); };
        auto nexts = [&](const char*& ref){ if (i+1 < argc) ref = argv[++i]; };
        if (a == "--w" || a == "-w") nexti(p.width);
        else if (a == "--h" || a == "-h") nexti(p.height);
        else if (a == "--fov" ) nextf(p.fov_y_deg);
        else if (a == "--theta") nextf(p.theta_deg);
        else if (a == "--mass") nextf(p.mass);
        else if (a == "--rcam") nextf(p.r_cam);
        else if (a == "--out") nexts(p.out_path);
        else if (a == "--samples" || a == "-s") nexti(p.samples);
        else if (a == "--spp-total") nexti(p.spp_total);
        else if (a == "--checkpoint-every") nexti(p.checkpoint_every);
        else if (a == "--spin") nextf(p.spin); // a/M
        else if (a == "--check-invariants") p.check_invariants = 1;
        else if (a == "--tile" ) nexti(p.tile_h);
        else if (a == "--frames") nexti(p.frames);
        else if (a == "--turns") nextf(p.turns);
        else if (a == "--phi0") nextf(p.phi0_deg);
        else if (a == "--outpat") nexts(p.out_pat);
        else if (a == "--mode" && i+1 < argc) p.mode = argv[++i];
    }
}

int main(int argc, char** argv) {
    RenderParams params;
    parse_args(argc, argv, params);

    int device_count = 0; cudaGetDeviceCount(&device_count);
    if (device_count <= 0) {
        std::fprintf(stderr, "No CUDA devices found.\n");
        return 1;
    }
    cudaSetDevice(0);

    auto deg2rad = [](float d){ return d * PI_F / 180.0f; };
    float theta_cam = deg2rad(params.theta_deg);
    float phi0 = deg2rad(params.phi0_deg);

    // Allocate reusable resources once
    const int tile_h = params.tile_h;
    const size_t tileN = (size_t)params.width * (size_t)tile_h;
    uchar3* d_tile[2] = {nullptr, nullptr};
    unsigned char* h_tile[2] = {nullptr, nullptr};
    cudaStream_t stream[2] = {nullptr, nullptr};
    float* d_inv[2] = {nullptr, nullptr};
    for (int i = 0; i < 2; ++i) {
        cudaMalloc(&d_tile[i], tileN * sizeof(uchar3));
        cudaHostAlloc((void**)&h_tile[i], tileN * 3, cudaHostAllocDefault);
        cudaStreamCreate(&stream[i]);
        if (params.check_invariants) {
            cudaMalloc(&d_inv[i], 3 * sizeof(float));
            cudaMemsetAsync(d_inv[i], 0, 3 * sizeof(float), stream[i]);
        }
    }

    int minGrid = 0, blockSize = 0;
    cudaOccupancyMaxPotentialBlockSize(&minGrid, &blockSize, renderKernel, 0, 0);
    if (blockSize <= 0) blockSize = 256;

    // Derive modes
    auto tolower_str = [](std::string s){ for (auto& c: s) c = (char)tolower(c); return s; };
    std::string mode = tolower_str(params.mode ? std::string(params.mode) : std::string("balanced"));
    int integrator_mode_host = 0; // 0=fixed, 1=adaptive
    int camera_mode_host = 0;     // 0=static,1=ZAMO
    if (mode == "accurate") { integrator_mode_host = 1; camera_mode_host = 1; params.check_invariants = 1; }
    else if (mode == "fast") { integrator_mode_host = 0; camera_mode_host = 0; }
    else { integrator_mode_host = 0; camera_mode_host = 0; }

    auto render_one = [&](float phi_cam, const char* path) -> bool {
        int total_spp = (params.spp_total > 0 ? params.spp_total : params.samples);
        int batch = std::max(1, params.samples);
        int passes = (total_spp + batch - 1) / batch;
        std::vector<float> accum((size_t)params.width * params.height * 3, 0.0f);
        int spp_done = 0;
        float inv_max_host[3] = {0,0,0};
        for (int p = 0; p < passes; ++p) {
            int pass_spp = std::min(batch, total_spp - spp_done);
            // Render all tiles for this pass and accumulate
            int tiles = (params.height + tile_h - 1) / tile_h;
            for (int t = 0; t < tiles + 1; ++t) {
                int buf = t % 2;
                int y0 = t * tile_h;
                int this_h = (y0 + tile_h <= params.height) ? tile_h : (params.height - y0);
                if (t < tiles) {
                    size_t thisN = (size_t)params.width * (size_t)tile_h;
                    dim3 block(blockSize);
                    dim3 grid((unsigned)((thisN + block.x - 1) / block.x));
                    // integrator_mode: 0=fixed RK4, 1=adaptive; camera_mode: 0=static, 1=ZAMO (host decides)
                    int integrator_mode = integrator_mode_host;
                    int camera_mode = camera_mode_host;
                    renderKernel<<<grid, block, 0, stream[buf]>>>(params.width, params.height, y0, tile_h, pass_spp,
                                                                 params.mass, params.spin, params.r_cam, theta_cam, phi_cam,
                                                                 params.fov_y_deg, params.check_invariants, integrator_mode, camera_mode,
                                                                 d_inv[buf], d_tile[buf]);
                    cudaMemcpyAsync(h_tile[buf], d_tile[buf], thisN * sizeof(uchar3), cudaMemcpyDeviceToHost, stream[buf]);
                }
                if (t > 0) {
                    int prev = (t - 1) % 2;
                    int prev_y0 = (t - 1) * tile_h;
                    int prev_h = (prev_y0 + tile_h <= params.height) ? tile_h : (params.height - prev_y0);
                    cudaStreamSynchronize(stream[prev]);
                    if (params.check_invariants) {
                        float inv3[3];
                        cudaMemcpy(inv3, d_inv[prev], 3 * sizeof(float), cudaMemcpyDeviceToHost);
                        inv_max_host[0] = std::max(inv_max_host[0], inv3[0]);
                        inv_max_host[1] = std::max(inv_max_host[1], inv3[1]);
                        inv_max_host[2] = std::max(inv_max_host[2], inv3[2]);
                        cudaMemsetAsync(d_inv[prev], 0, 3 * sizeof(float), stream[prev]);
                    }
                    // Accumulate this tile
                    for (int row = 0; row < prev_h; ++row) {
                        int y = prev_y0 + row;
                        size_t base_img = (size_t)y * params.width * 3;
                        size_t base_tile = (size_t)row * params.width * 3;
                        for (int x = 0; x < params.width; ++x) {
                            unsigned char r = h_tile[prev][base_tile + 3*x + 0];
                            unsigned char g = h_tile[prev][base_tile + 3*x + 1];
                            unsigned char b = h_tile[prev][base_tile + 3*x + 2];
                            accum[base_img + 3*x + 0] += (float)r / 255.0f;
                            accum[base_img + 3*x + 1] += (float)g / 255.0f;
                            accum[base_img + 3*x + 2] += (float)b / 255.0f;
                        }
                    }
                }
            }
            spp_done += pass_spp;
            bool do_ckpt = (p == passes - 1) || (params.checkpoint_every > 0 && ((p+1) % params.checkpoint_every == 0));
            if (do_ckpt) {
                std::vector<unsigned char> out((size_t)params.width * params.height * 3);
                float inv = 1.0f / (float)spp_done;
                for (size_t i = 0, n = out.size()/3; i < n; ++i) {
                    float r = fminf(fmaxf(accum[3*i + 0] * inv, 0.0f), 1.0f);
                    float g = fminf(fmaxf(accum[3*i + 1] * inv, 0.0f), 1.0f);
                    float b = fminf(fmaxf(accum[3*i + 2] * inv, 0.0f), 1.0f);
                    out[3*i + 0] = (unsigned char)(255.0f * r);
                    out[3*i + 1] = (unsigned char)(255.0f * g);
                    out[3*i + 2] = (unsigned char)(255.0f * b);
                }
                if (!write_ppm(path, params.width, params.height, out)) {
                    std::fprintf(stderr, "Failed to write %s\n", path); return false;
                }
                if (params.check_invariants) {
                    std::printf("[checkpoint] %s @ %d spp  H|max=%.3e  dE/E|max=%.3e  dL/L|max=%.3e\n",
                                path, spp_done, inv_max_host[0], inv_max_host[1], inv_max_host[2]);
                } else {
                    std::printf("[checkpoint] %s @ %d spp\n", path, spp_done);
                }
            }
        }
        return true;
    };

    if (params.frames <= 1) {
        if (!render_one(phi0, params.out_path)) return 1;
        std::printf("Wrote %s\n", params.out_path);
    } else {
        for (int i = 0; i < params.frames; ++i) {
            float t = (float)i / (float)params.frames;
            float phi = phi0 + 2.0f * PI_F * params.turns * t;
            char path[512];
            snprintf(path, sizeof(path), params.out_pat, i);
            if (!render_one(phi, path)) return 1;
        }
        std::printf("Wrote %d frames matching pattern %s\n", params.frames, params.out_pat);
    }

    for (int i = 0; i < 2; ++i) {
        cudaFree(d_tile[i]);
        cudaFreeHost(h_tile[i]);
        cudaStreamDestroy(stream[i]);
        if (d_inv[i]) cudaFree(d_inv[i]);
    }
    return 0;
}
