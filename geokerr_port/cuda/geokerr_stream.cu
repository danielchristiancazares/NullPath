/**
 * CUDA per-λ stream sampler for Kerr photon geodesics (Mino time λ).
 *
 * Generates u(λ), mu(λ), t(λ), φ(λ), and affine(λ) over a uniform λ-grid
 * using RK2 (midpoint) integration of the correct equations of motion:
 *  - dr/dλ = ±sqrt(R(r))
 *  - dμ/dλ = ±sqrt((1-μ^2) Q - μ^2 Lz^2)
 *  - dφ/dλ, dt/dλ from standard Kerr photon formulas (E=1)
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>

#ifndef GKSTREAM_TINY
#define GKSTREAM_TINY 1e-15
#endif

__device__ inline double delta_fn(double r, double a) {
    return r*r - 2.0*r + a*a;
}

// Photon radial potential R(r)
__device__ inline double R_of_r(double r, double a, double Lz, double Q2) {
    // E=1
    double Delta = delta_fn(r, a);
    double term1 = (r*r + a*a) - a*Lz; term1 *= term1;
    double term2 = Delta * (Q2 + (Lz - a)*(Lz - a));
    return term1 - term2;
}

// (dμ/dλ)^2 for photons in Mino time
__device__ inline double Theta_mu(double mu, double a, double Lz, double Q2) {
    (void)a; // a^2(1-E^2) term vanishes for E=1
    double one_minus_mu2 = fmax(0.0, 1.0 - mu*mu);
    return one_minus_mu2 * Q2 - (mu*mu) * (Lz*Lz);
}

// dφ/dλ (E=1), unwrapped
__device__ inline double dphi_dlam(double r, double mu, double a, double Lz) {
    double Delta = fmax(GKSTREAM_TINY, delta_fn(r, a));
    double sin2 = fmax(GKSTREAM_TINY, 1.0 - mu*mu);
    return a * (((r*r + a*a) - a*Lz) / Delta) + Lz / sin2 - a;
}

// dt/dλ (E=1)
__device__ inline double dt_dlam(double r, double mu, double a, double Lz) {
    double Delta = fmax(GKSTREAM_TINY, delta_fn(r, a));
    return ((r*r + a*a) * ((r*r + a*a) - a*Lz)) / Delta + a * (Lz - a * (1.0 - mu*mu));
}

__device__ inline double clamp(double x, double lo, double hi) {
    return fmin(hi, fmax(lo, x));
}

// Kernel: single-thread RK4 stream integration with sub-steps over S samples
__global__ void geokerr_stream_kernel(
    double u0, double uf, double mu0, double a, double L, double Q2,
    double lam0, double dlam, int S,
    double* u, double* mu, double* tt, double* phi, double* affine)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;

    // Initialize state
    double r = (u0 > GKSTREAM_TINY) ? 1.0 / u0 : 1.0 / GKSTREAM_TINY;
    // Keep outside the outer horizon
    double rh = 1.0 + sqrt(fmax(0.0, 1.0 - a*a));
    if (r < rh + 1e-6) r = rh + 1e-6;
    double m = clamp(mu0, -1.0, 1.0);
    double t = 0.0;
    double p = 0.0;

    // Initial direction: Fortran SU is sign in u-space; dr/dλ sign is opposite of SU
    int su = (uf > u0) ? +1 : (uf < u0 ? -1 : +1);
    int sr = -su;
    int sm = (mu0 <= 0.0) ? +1 : -1;

    // Sub-stepping configuration
    const int NSUB = 64; // Increased from 16 to 64 for better accuracy
    const double h = (NSUB > 0) ? (dlam / (double)NSUB) : dlam;

    for (int i = 0; i < S; i++) {
        double lam = lam0 + (double)i * dlam;
        // Write outputs
        affine[i] = lam;
        // Force exact seeding for first sample to match inputs
        if (i == 0) {
            u[i] = u0;
            mu[i] = mu0;
        } else {
            u[i] = (r > GKSTREAM_TINY) ? 1.0 / r : 1.0 / GKSTREAM_TINY;
            mu[i] = m;
        }
        tt[i] = t;
        phi[i] = p;

        if (i == S - 1) break;
        // Advance by dlam using NSUB RK4 sub-steps with sign/turning-point handling
        for (int s = 0; s < NSUB; ++s) {
            // k1
            double R1 = R_of_r(r, a, L, Q2);
            if (R1 < 0.0) { R1 = 0.0; sr = -sr; }
            double T1 = Theta_mu(m, a, L, Q2);
            if (T1 < 0.0) { T1 = 0.0; sm = -sm; }
            double k1_r = sr * sqrt(R1);
            double k1_m = sm * sqrt(T1);
            double k1_p = dphi_dlam(r, m, a, L);
            double k1_t = dt_dlam(r, m, a, L);

            // k2
            double r2 = r + 0.5 * h * k1_r;
            double m2 = clamp(m + 0.5 * h * k1_m, -1.0, 1.0);
            double R2 = R_of_r(r2, a, L, Q2);
            if (R2 < 0.0) { R2 = 0.0; sr = -sr; r2 = r; }
            double T2 = Theta_mu(m2, a, L, Q2);
            if (T2 < 0.0) { T2 = 0.0; sm = -sm; m2 = m; }
            double k2_r = sr * sqrt(R2);
            double k2_m = sm * sqrt(T2);
            double k2_p = dphi_dlam(r2, m2, a, L);
            double k2_t = dt_dlam(r2, m2, a, L);

            // k3
            double r3 = r + 0.5 * h * k2_r;
            double m3 = clamp(m + 0.5 * h * k2_m, -1.0, 1.0);
            double R3 = R_of_r(r3, a, L, Q2);
            if (R3 < 0.0) { R3 = 0.0; sr = -sr; r3 = r; }
            double T3 = Theta_mu(m3, a, L, Q2);
            if (T3 < 0.0) { T3 = 0.0; sm = -sm; m3 = m; }
            double k3_r = sr * sqrt(R3);
            double k3_m = sm * sqrt(T3);
            double k3_p = dphi_dlam(r3, m3, a, L);
            double k3_t = dt_dlam(r3, m3, a, L);

            // k4
            double r4 = r + h * k3_r;
            double m4 = clamp(m + h * k3_m, -1.0, 1.0);
            double R4 = R_of_r(r4, a, L, Q2);
            if (R4 < 0.0) { R4 = 0.0; sr = -sr; r4 = r; }
            double T4 = Theta_mu(m4, a, L, Q2);
            if (T4 < 0.0) { T4 = 0.0; sm = -sm; m4 = m; }
            double k4_r = sr * sqrt(R4);
            double k4_m = sm * sqrt(T4);
            double k4_p = dphi_dlam(r4, m4, a, L);
            double k4_t = dt_dlam(r4, m4, a, L);

            // RK4 update
            r += (h/6.0) * (k1_r + 2.0*k2_r + 2.0*k3_r + k4_r);
            m = clamp(m + (h/6.0) * (k1_m + 2.0*k2_m + 2.0*k3_m + k4_m), -1.0, 1.0);
            p += (h/6.0) * (k1_p + 2.0*k2_p + 2.0*k3_p + k4_p);
            t += (h/6.0) * (k1_t + 2.0*k2_t + 2.0*k3_t + k4_t);

            // Enforce horizon clamp
            r = fmax(rh + 1e-12, r);
        }
    }
}

extern "C" int geokerr_stream_cuda(
    double u0, double uf, double mu0, double a, double L, double Q2,
    double lam0, double dlam, int S,
    double* h_u, double* h_mu, double* h_t, double* h_phi, double* h_affine)
{
    if (S <= 0 || !h_u || !h_mu || !h_t || !h_phi || !h_affine) return -1;

    // Allocate device buffers
    double *d_u = nullptr, *d_mu = nullptr, *d_t = nullptr, *d_phi = nullptr, *d_aff = nullptr;
    size_t bytes = (size_t)S * sizeof(double);
    cudaError_t err;
    if ((err = cudaMalloc(&d_u, bytes)) != cudaSuccess) return -2;
    if ((err = cudaMalloc(&d_mu, bytes)) != cudaSuccess) { cudaFree(d_u); return -2; }
    if ((err = cudaMalloc(&d_t, bytes)) != cudaSuccess) { cudaFree(d_u); cudaFree(d_mu); return -2; }
    if ((err = cudaMalloc(&d_phi, bytes)) != cudaSuccess) { cudaFree(d_u); cudaFree(d_mu); cudaFree(d_t); return -2; }
    if ((err = cudaMalloc(&d_aff, bytes)) != cudaSuccess) { cudaFree(d_u); cudaFree(d_mu); cudaFree(d_t); cudaFree(d_phi); return -2; }

    // Launch single-thread kernel (one geodesic stream)
    geokerr_stream_kernel<<<1,1>>>(u0, uf, mu0, a, L, Q2, lam0, dlam, S, d_u, d_mu, d_t, d_phi, d_aff);
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        cudaFree(d_u); cudaFree(d_mu); cudaFree(d_t); cudaFree(d_phi); cudaFree(d_aff);
        return -3;
    }

    // Copy back
    cudaMemcpy(h_u, d_u, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mu, d_mu, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_t, d_t, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_phi, d_phi, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_affine, d_aff, bytes, cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_u); cudaFree(d_mu); cudaFree(d_t); cudaFree(d_phi); cudaFree(d_aff);
    return 0;
}
