// schwarzschild_blackhole.cu
// NullPath — Schwarzschild Null‑Geodesic Tracer (CUDA)
// Compile: nvcc -O3 -arch=sm_89 schwarzschild_blackhole.cu -o bin/nullpath

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <chrono>
#include <cstring>
#include <cstdlib>

// Centralized constants and shared types
#include "include/bh_common.h"
#include "include/bh_types.h"
#include "include/bh_schwarzschild.h"

// (Types moved to include/bh_types.h)


// Analytic turning radius r_min(b) for Schwarzschild equatorial photons
// Solves r^3 - b^2 r + b^2 Rs = 0 using trigonometric closed form.
// Returns true and sets rmin when b > b_crit; returns false if no turning point (b <= b_crit).
__device__ __host__ inline bool schwarzschild_rmin_analytic(float Rs, float b, float* rmin_out) {
    if (Rs <= 0.0f || b <= 0.0f || rmin_out == nullptr) return false;
    const double beta = (double)b / (double)Rs; // dimensionless b/Rs
    const double beta_crit = 0.5 * 3.0 * sqrt(3.0); // (3√3/2)
    if (beta <= beta_crit) return false; // no turning point outside horizon
    const double p = -beta * beta;
    const double A = 2.0 * sqrt(-p / 3.0);         // = 2 beta / sqrt(3)
    double cosarg = -(3.0 * sqrt(3.0)) / (2.0 * beta); // = -beta_crit/beta ∈ [-1, 0)
    if (cosarg < -1.0) cosarg = -1.0; if (cosarg > 1.0) cosarg = 1.0;
    const double ang = (1.0 / 3.0) * acos(cosarg);
    const double u = A * cos(ang); // largest real root in u = r/Rs
    *rmin_out = (float)(u * (double)Rs);
    return isfinite(*rmin_out) && *rmin_out > Rs;
}

// Schwarzschild geodesic integration lives in include/bh_schwarzschild.h

// (GeodesicLookupTable moved to include/bh_types.h)

// CUDA kernel to precompute geodesic deflection data
__global__ 
void precomputeGeodesicsKernel(float mass, float r_start, float* impact_parameters, 
                              float* deflection_angles, float* closest_approaches,
                              int* ray_status, int num_rays) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_rays) return;
    
    SchwarzschildMetric metric(mass);
    float b = impact_parameters[idx]; // Impact parameter
    
    // Initialize geodesic for light ray from infinity
    GeodesicState state;
    state.r = r_start;  // Start from "infinity"
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
        if (state.r > r_start + 1.0f && state.pr > 0) {
            escaped = true;
            break;
        }
    }
    
    // Store results (with turning-point refinement for escaped rays)
    deflection_angles[idx] = total_phi_change;
    float rmin_out = min_r;
    if (escaped) {
        float r_an = 0.0f;
        if (schwarzschild_rmin_analytic(metric.Rs, b, &r_an)) {
            rmin_out = r_an; // exact turning radius (GeoKerr convention)
        }
    }
    closest_approaches[idx] = rmin_out;
    
    if (escaped || (state.r > r_start + 1.0f && state.pr > 0)) {
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
    if (num_rays < 2) {
        std::fprintf(stderr, "num_rays must be >= 2.\n");
        return nullptr;
    }
    // Reference values for comparison and GeoKerr alignment
    const float Rs = SCHWARZSCHILD_MULTIPLIER * mass;
    const float bcrit_Rs = bh_bcrit_from_Rs(Rs);
    const float bcrit_M  = bh_bcrit_from_M(mass);
    printf("[ref] r_ph = %.6f (%.3f Rs), b_crit = %.6f (%.6f M)\n",
           bh_r_ph(Rs), BH_PHOTON_SPHERE_MULTIPLIER, bcrit_Rs, bcrit_M);
    
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
    float b_min = 0.30f * bcrit_Rs;   // Focus resolution around the critical region
    float b_max = 20.0f * bcrit_Rs;   // Still spans far field
    if (const char* bmin_env = std::getenv("BH_B_MIN_RS")) {
        float v = strtof(bmin_env, nullptr); if (v > 0.f) b_min = v * Rs;
    }
    if (const char* bmax_env = std::getenv("BH_B_MAX_RS")) {
        float v = strtof(bmax_env, nullptr); if (v > 0.f) b_max = v * Rs;
    }
    if (b_min >= b_max) { b_min = Rs; b_max = 50.0f * Rs; }
    
    for (int i = 0; i < num_rays; i++) {
        float t = (float)i / (num_rays - 1);
        // Logarithmic spacing for better resolution near critical values
        h_impact_params[i] = b_min * powf(b_max / b_min, t);
    }

    float r_start = 1000.0f;
    float min_start = 1.2f * b_max;
    if (min_start > r_start) r_start = min_start;
    
    // Allocate device memory
    float *d_impact_params = nullptr, *d_deflection_angles = nullptr, *d_closest_approaches = nullptr;
    int *d_ray_status = nullptr;
    auto cleanup_device = [&]() {
        if (d_impact_params) cudaFree(d_impact_params);
        if (d_deflection_angles) cudaFree(d_deflection_angles);
        if (d_closest_approaches) cudaFree(d_closest_approaches);
        if (d_ray_status) cudaFree(d_ray_status);
    };

    cudaError_t err;
    err = cudaMalloc(&d_impact_params, num_rays * sizeof(float));
    if (err) { fprintf(stderr, "cudaMalloc d_impact_params failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    err = cudaMalloc(&d_deflection_angles, num_rays * sizeof(float));
    if (err) { fprintf(stderr, "cudaMalloc d_deflection_angles failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    err = cudaMalloc(&d_closest_approaches, num_rays * sizeof(float));
    if (err) { fprintf(stderr, "cudaMalloc d_closest_approaches failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    err = cudaMalloc(&d_ray_status, num_rays * sizeof(int));
    if (err) { fprintf(stderr, "cudaMalloc d_ray_status failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }

    // Copy input data to device
    err = cudaMemcpy(d_impact_params, h_impact_params.data(), num_rays * sizeof(float), cudaMemcpyHostToDevice);
    if (err) { fprintf(stderr, "cudaMemcpy H2D impact_params failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    
    // Launch kernel
    int block_size = 256;
    int num_blocks = (num_rays + block_size - 1) / block_size;
    
    precomputeGeodesicsKernel<<<num_blocks, block_size>>>(
        mass, r_start, d_impact_params, d_deflection_angles, d_closest_approaches, 
        d_ray_status, num_rays);
    
    // Check for launch errors and synchronize
    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(launch_err));
        cleanup_device(); delete table; return nullptr;
    }

    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "Kernel execution failed: %s\n", cudaGetErrorString(err));
        cleanup_device(); delete table; return nullptr;
    }
    
    // Copy results back to host
    err = cudaMemcpy(h_deflection_angles.data(), d_deflection_angles, num_rays * sizeof(float), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H deflection_angles failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    err = cudaMemcpy(h_closest_approaches.data(), d_closest_approaches, num_rays * sizeof(float), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H closest_approaches failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    err = cudaMemcpy(h_ray_status.data(), d_ray_status, num_rays * sizeof(int), cudaMemcpyDeviceToHost);
    if (err) { fprintf(stderr, "cudaMemcpy D2H ray_status failed: %s\n", cudaGetErrorString(err)); cleanup_device(); delete table; return nullptr; }
    
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
    cleanup_device();
    
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
