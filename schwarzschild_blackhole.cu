// schwarzschild_blackhole.cu
// NullPath — Schwarzschild Null‑Geodesic Tracer (CUDA)
// Compile: nvcc -O3 -arch=sm_89 schwarzschild_blackhole.cu -o bin/nullpath

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <chrono>
#include <cstring>

// Centralized constants and shared types
#include "include/bh_common.h"
#include "include/bh_types.h"

// (Types moved to include/bh_types.h)

// Geodesic equations of motion in Schwarzschild spacetime
__device__
void computeGeodesicDerivatives(const GeodesicState& state, const SchwarzschildMetric& metric, 
                               GeodesicState& derivatives) {
    float r = state.r;
    float theta = state.theta;
    float Rs = metric.Rs;
    
    // Avoid singularities
    if (r <= Rs + 0.001f || r <= 0.001f) {
        // Inside event horizon or at singularity - trajectory ends
        derivatives = GeodesicState();
        return;
    }
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float r2 = r * r;
    float r3 = r2 * r;
    float A = 1.0f - Rs / r; // g_rr^{-1} = g^{rr}

    // Hamiltonian form: dx^μ/dλ = g^{μν} p_ν
    // dt/dλ = g^{tt} p_t = (-1/A) * p_t
    derivatives.t = (-1.0f / A) * state.pt;

    // dr/dλ = g^{rr} p_r = A * p_r
    derivatives.r = A * state.pr;

    // dθ/dλ = g^{θθ} p_θ = (1/r^2) * p_θ
    derivatives.theta = state.ptheta / r2;

    // dφ/dλ = g^{φφ} p_φ = (1/(r^2 sin^2θ)) * p_φ
    derivatives.phi = state.pphi / (r2 * fmaxf(1e-12f, sin_theta * sin_theta));

    // 4-momentum derivatives: dp_μ/dλ = -½ ∂_μ g^{αβ} p_α p_β
    derivatives.pt = 0.0f; // t is cyclic (metric independent of t)

    // r-derivatives of inverse metric components
    float sin2 = fmaxf(1e-12f, sin_theta * sin_theta);
    float dgtt_inv_dr = Rs / (r2 * A * A);               // ∂_r g^{tt}
    float dgrr_inv_dr = Rs / (r2);                       // ∂_r g^{rr}
    float dgthth_inv_dr = -2.0f / (r3);                  // ∂_r g^{θθ}
    float dgphph_inv_dr = -2.0f / (r3 * sin2);           // ∂_r g^{φφ}

    derivatives.pr = -0.5f * (
        dgtt_inv_dr   * state.pt     * state.pt +
        dgrr_inv_dr   * state.pr     * state.pr +
        dgthth_inv_dr * state.ptheta * state.ptheta +
        dgphph_inv_dr * state.pphi   * state.pphi
    );

    // θ-momentum: only g^{φφ} depends on θ
    // Protect sin^3(theta) in denominator with sign-preserving clamp
    float sin_theta_safe = (fabsf(sin_theta) < 1e-6f) ? (sin_theta >= 0.0f ? 1e-6f : -1e-6f) : sin_theta;
    derivatives.ptheta = (cos_theta / (r2 * sin2 * sin_theta_safe)) * state.pphi * state.pphi;

    // φ is cyclic
    derivatives.pphi = 0.0f;
}

// 4th order Runge-Kutta integration for geodesics
__device__
bool integrateGeodesic(GeodesicState& state, const SchwarzschildMetric& metric, 
                      float step_size) {
    GeodesicState k1, k2, k3, k4;
    GeodesicState temp_state;
    
    // k1 = f(state)
    computeGeodesicDerivatives(state, metric, k1);
    
    // k2 = f(state + 0.5*h*k1)
    temp_state.t = state.t + 0.5f * step_size * k1.t;
    temp_state.r = state.r + 0.5f * step_size * k1.r;
    temp_state.theta = state.theta + 0.5f * step_size * k1.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k1.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k1.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k1.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k1.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k1.pphi;
    
    computeGeodesicDerivatives(temp_state, metric, k2);
    
    // k3 = f(state + 0.5*h*k2)
    temp_state.t = state.t + 0.5f * step_size * k2.t;
    temp_state.r = state.r + 0.5f * step_size * k2.r;
    temp_state.theta = state.theta + 0.5f * step_size * k2.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k2.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k2.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k2.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k2.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k2.pphi;
    
    computeGeodesicDerivatives(temp_state, metric, k3);
    
    // k4 = f(state + h*k3)
    temp_state.t = state.t + step_size * k3.t;
    temp_state.r = state.r + step_size * k3.r;
    temp_state.theta = state.theta + step_size * k3.theta;
    temp_state.phi = state.phi + step_size * k3.phi;
    temp_state.pt = state.pt + step_size * k3.pt;
    temp_state.pr = state.pr + step_size * k3.pr;
    temp_state.ptheta = state.ptheta + step_size * k3.ptheta;
    temp_state.pphi = state.pphi + step_size * k3.pphi;
    
    computeGeodesicDerivatives(temp_state, metric, k4);
    
    // Update state: state = state + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    float h_6 = step_size / 6.0f;
    state.t += h_6 * (k1.t + 2.0f*k2.t + 2.0f*k3.t + k4.t);
    state.r += h_6 * (k1.r + 2.0f*k2.r + 2.0f*k3.r + k4.r);
    state.theta += h_6 * (k1.theta + 2.0f*k2.theta + 2.0f*k3.theta + k4.theta);
    state.phi += h_6 * (k1.phi + 2.0f*k2.phi + 2.0f*k3.phi + k4.phi);
    state.pt += h_6 * (k1.pt + 2.0f*k2.pt + 2.0f*k3.pt + k4.pt);
    state.pr += h_6 * (k1.pr + 2.0f*k2.pr + 2.0f*k3.pr + k4.pr);
    state.ptheta += h_6 * (k1.ptheta + 2.0f*k2.ptheta + 2.0f*k3.ptheta + k4.ptheta);
    state.pphi += h_6 * (k1.pphi + 2.0f*k2.pphi + 2.0f*k3.pphi + k4.pphi);
    
    // Check for termination conditions (allow r to grow beyond start to detect escape)
    return (state.r > metric.Rs + 0.001f &&
            state.theta > 0.001f && state.theta < PI_F - 0.001f);
}

// (GeodesicLookupTable moved to include/bh_types.h)

// CUDA kernel to precompute geodesic deflection data
__global__ 
void precomputeGeodesicsKernel(float mass, float* impact_parameters, 
                              float* deflection_angles, float* closest_approaches,
                              int* ray_status, int num_rays) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_rays) return;
    
    SchwarzschildMetric metric(mass);
    float b = impact_parameters[idx]; // Impact parameter
    
    // Initialize geodesic for light ray from infinity
    GeodesicState state;
    state.r = 1000.0f;  // Start from "infinity"
    state.theta = PI_F / 2.0f;  // Equatorial plane
    state.phi = 0.0f;
    state.t = 0.0f;
    
    // Initial 4-momentum for photon (null geodesic)
    float energy = 1.0f;            // Photon energy at infinity (E)
    state.pt = -energy;             // p_t = -E (conserved)
    state.pphi = b * energy;        // p_φ = L = bE (conserved)
    state.ptheta = 0.0f;            // Equatorial plane

    // Radial momentum from null Hamiltonian: H = ½ g^{μν} p_μ p_ν = 0
    float r = state.r;
    float A = 1.0f - metric.Rs / r;
    float L = state.pphi;
    float pr_squared = (energy * energy) / (A * A) - (L * L) / (A * r * r);
    state.pr = -sqrtf(fmaxf(0.0f, pr_squared)); // Negative: ingoing from infinity
    
    float min_r = state.r;
    float total_phi_change = 0.0f;
    float prev_phi = state.phi;
    bool escaped = false;
    float prev_r = state.r;
    float prev_pr = state.pr;
    
    // Integrate geodesic with radius-aware step scaling (larger steps far away)
    for (int step = 0; step < MAX_INTEGRATION_STEPS; step++) {
        // Aggressive far-field stepping: up to 100x base step when r >> Rs
        float scale = fminf(100.0f, fmaxf(1.0f, state.r / (0.5f * metric.Rs)));
        float h = INTEGRATION_STEP_SIZE * scale;
        if (!integrateGeodesic(state, metric, h)) {
            break; // final classification done after loop
        }
        
        // Track minimum approach (refined at turning point)
        min_r = fminf(min_r, state.r);
        if (prev_pr < 0.0f && state.pr > 0.0f) {
            float denom = (state.pr - prev_pr);
            if (fabsf(denom) > 1e-12f) {
                float alpha = (-prev_pr) / denom; // fraction from prev to curr where p_r=0
                alpha = fminf(fmaxf(alpha, 0.0f), 1.0f);
                float r_turn = prev_r + alpha * (state.r - prev_r);
                min_r = fminf(min_r, r_turn);
            }
        }
        prev_r = state.r;
        prev_pr = state.pr;
        
        // Track total deflection angle
        float phi_diff = state.phi - prev_phi;
        // Handle phi wraparound
        if (phi_diff > PI_F) phi_diff -= 2.0f * PI_F;
        if (phi_diff < -PI_F) phi_diff += 2.0f * PI_F;
        total_phi_change += phi_diff;
        prev_phi = state.phi;
        
        // Check if ray escaped to infinity
        if (state.r > 1001.0f && state.pr > 0) {
            escaped = true;
            break;
        }
    }
    
    // Store results (with turning-point refinement for escaped rays)
    deflection_angles[idx] = total_phi_change;
    float rmin_out = min_r;
    if (escaped && min_r > metric.Rs + 1e-4f) {
        // For Schwarzschild equatorial photons, turning radius r satisfies:
        // r^3 - b^2 r + b^2 Rs = 0 (with b = impact parameter).
        // Do a couple of Newton steps starting from min_r.
        float r_guess = fmaxf(min_r, metric.Rs + 1e-3f);
        float b2 = b * b;
        #pragma unroll 2
        for (int it = 0; it < 2; ++it) {
            float f = r_guess*r_guess*r_guess - b2 * r_guess + b2 * metric.Rs;
            float fp = 3.0f * r_guess * r_guess - b2;
            if (fabsf(fp) < 1e-8f) break;
            float step = f / fp;
            r_guess = fmaxf(metric.Rs + 1e-3f, r_guess - step);
        }
        // Keep the smaller of integrated min and refined value
        if (isfinite(r_guess)) rmin_out = fminf(min_r, r_guess);
    }
    closest_approaches[idx] = rmin_out;
    
    if (escaped || (state.r > 1001.0f && state.pr > 0)) {
        ray_status[idx] = 0; // Escaped to infinity (definitive)
    } else if (state.r <= metric.Rs + 0.01f) {
        ray_status[idx] = 1; // Captured by black hole
    } else if (fabsf(min_r - PHOTON_SPHERE_MULTIPLIER * metric.Rs) < 0.01f * metric.Rs) {
        ray_status[idx] = 2; // Near photon sphere
    } else {
        ray_status[idx] = 3; // Inconclusive within step budget
    }
}

// Host function to precompute geodesic lookup table
GeodesicLookupTable* precomputeGeodesics(float mass, int num_rays) {
    printf("Precomputing %d geodesics for black hole mass %.2f...\n", num_rays, mass);
    
    // Allocate host memory
    GeodesicLookupTable* table = new GeodesicLookupTable();
    table->mass = mass;
    table->num_entries = num_rays;
    
    std::vector<float> h_impact_params(num_rays);
    std::vector<float> h_deflection_angles(num_rays);
    std::vector<float> h_closest_approaches(num_rays);
    std::vector<float> h_bending_angles(num_rays);
    std::vector<int> h_ray_status(num_rays);
    
    // Generate impact parameter range
    float Rs = SCHWARZSCHILD_MULTIPLIER * mass;
    float b_min = Rs; // Start just outside Schwarzschild radius
    float b_max = 50.0f * Rs; // Extended range for strong lensing
    
    for (int i = 0; i < num_rays; i++) {
        float t = (float)i / (num_rays - 1);
        // Logarithmic spacing for better resolution near critical values
        h_impact_params[i] = b_min * powf(b_max / b_min, t);
    }
    
    // Allocate device memory
    float *d_impact_params, *d_deflection_angles, *d_closest_approaches;
    int *d_ray_status;
    
    cudaError_t err;
    err = cudaMalloc(&d_impact_params, num_rays * sizeof(float)); if (err) { fprintf(stderr, "cudaMalloc d_impact_params failed: %s\n", cudaGetErrorString(err)); return nullptr; }
    err = cudaMalloc(&d_deflection_angles, num_rays * sizeof(float)); if (err) { fprintf(stderr, "cudaMalloc d_deflection_angles failed: %s\n", cudaGetErrorString(err)); cudaFree(d_impact_params); return nullptr; }
    err = cudaMalloc(&d_closest_approaches, num_rays * sizeof(float)); if (err) { fprintf(stderr, "cudaMalloc d_closest_approaches failed: %s\n", cudaGetErrorString(err)); cudaFree(d_impact_params); cudaFree(d_deflection_angles); return nullptr; }
    err = cudaMalloc(&d_ray_status, num_rays * sizeof(int)); if (err) { fprintf(stderr, "cudaMalloc d_ray_status failed: %s\n", cudaGetErrorString(err)); cudaFree(d_impact_params); cudaFree(d_deflection_angles); cudaFree(d_closest_approaches); return nullptr; }
    
    // Copy input data to device
    err = cudaMemcpy(d_impact_params, h_impact_params.data(), num_rays * sizeof(float), cudaMemcpyHostToDevice);
    if (err) { fprintf(stderr, "cudaMemcpy H2D impact_params failed: %s\n", cudaGetErrorString(err)); }
    
    // Launch kernel
    int block_size = 256;
    int num_blocks = (num_rays + block_size - 1) / block_size;
    
    precomputeGeodesicsKernel<<<num_blocks, block_size>>>(
        mass, d_impact_params, d_deflection_angles, d_closest_approaches, 
        d_ray_status, num_rays);
    
    // Check for launch errors and synchronize
    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(launch_err));
    }
    
    cudaDeviceSynchronize();
    cudaError_t sync_err = cudaGetLastError();
    if (sync_err != cudaSuccess) {
        fprintf(stderr, "Kernel execution failed: %s\n", cudaGetErrorString(sync_err));
    }
    
    // Copy results back to host
    err = cudaMemcpy(h_deflection_angles.data(), d_deflection_angles, num_rays * sizeof(float), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H deflection_angles failed: %s\n", cudaGetErrorString(err)); }
    err = cudaMemcpy(h_closest_approaches.data(), d_closest_approaches, num_rays * sizeof(float), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H closest_approaches failed: %s\n", cudaGetErrorString(err)); }
    err = cudaMemcpy(h_ray_status.data(), d_ray_status, num_rays * sizeof(int), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H ray_status failed: %s\n", cudaGetErrorString(err)); }
    
    // Allocate and copy final results
    table->impact_parameters = new float[num_rays];
    table->deflection_angles = new float[num_rays];
    table->closest_approaches = new float[num_rays]; 
    table->bending_angles = new float[num_rays];
    table->ray_status = new int[num_rays];
    
    memcpy(table->impact_parameters, h_impact_params.data(), num_rays * sizeof(float));
    memcpy(table->deflection_angles, h_deflection_angles.data(), num_rays * sizeof(float));
    memcpy(table->closest_approaches, h_closest_approaches.data(), num_rays * sizeof(float));
    // Compute bending alpha = total_phi_change - pi for escaped rays; NaN otherwise
    for (int i = 0; i < num_rays; ++i) {
        if (h_ray_status[i] == 0) {
            h_bending_angles[i] = h_deflection_angles[i] - PI_F;
        } else {
            h_bending_angles[i] = NAN;
        }
    }
    memcpy(table->bending_angles, h_bending_angles.data(), num_rays * sizeof(float));
    memcpy(table->ray_status, h_ray_status.data(), num_rays * sizeof(int));
    
    // Clean up device memory
    cudaFree(d_impact_params);
    cudaFree(d_deflection_angles);
    cudaFree(d_closest_approaches);
    cudaFree(d_ray_status);
    
    printf("Geodesic precomputation complete!\n");
    return table;
}

// Simple test program
#ifndef BLACKHOLE_NO_MAIN
int main() {
    printf("NullPath — CUDA null-geodesic tracer (Schwarzschild)\n");
    printf("====================================================\n");
    
    // Initialize CUDA
    cudaSetDevice(0);
    
    // Print device info
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Using GPU: %s\n", prop.name);
    printf("CUDA Cores (approx): ~%d\n", prop.multiProcessorCount * 128);
    printf("Global Memory: %.1f GB\n", prop.totalGlobalMem / (1024.0*1024.0*1024.0));
    
    // Test parameters
    float black_hole_mass = 10.0f; // Solar masses
    int num_geodesics = 100000;    // default resolution
    
    // Precompute geodesics
    auto start = std::chrono::high_resolution_clock::now();
    GeodesicLookupTable* table = precomputeGeodesics(black_hole_mass, num_geodesics);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    printf("Precomputation took %ld ms\n", duration.count());
    
    // Analysis results
    int captured = 0, escaped = 0, photon_sphere = 0;
    for (int i = 0; i < num_geodesics; i++) {
        switch (table->ray_status[i]) {
            case 0: escaped++; break;
            case 1: captured++; break;
            case 2: photon_sphere++; break;
        }
    }
    
    printf("\nResults for mass %.2f solar masses:\n", black_hole_mass);
    printf("Escaped rays: %d (%.1f%%)\n", escaped, 100.0f * escaped / num_geodesics);
    printf("Captured rays: %d (%.1f%%)\n", captured, 100.0f * captured / num_geodesics);
    printf("Photon sphere: %d (%.1f%%)\n", photon_sphere, 100.0f * photon_sphere / num_geodesics);
    
    // Clean up
    delete[] table->impact_parameters;
    delete[] table->deflection_angles;
    delete[] table->closest_approaches;
    delete[] table->bending_angles;
    delete[] table->ray_status;
    delete table;
    
    return 0;
}
#endif // BLACKHOLE_NO_MAIN
