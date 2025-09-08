// Shared physics types (header-only, inlined for device/host use).
// Intent: centralize structs used by both precompute and render paths.

#ifndef BH_TYPES_H
#define BH_TYPES_H

#include "bh_common.h"

// Schwarzschild metric helpers (geometric units, G=c=1)
struct SchwarzschildMetric {
    float mass;  // M
    float Rs;    // 2M

    __host__ __device__
    explicit SchwarzschildMetric(float m) : mass(m), Rs(BH_SCHWARZSCHILD_MULTIPLIER * m) {}

    __host__ __device__ inline float g_tt(float r) const { return -(1.0f - Rs / r); }
    __host__ __device__ inline float g_rr(float r) const { return  1.0f / (1.0f - Rs / r); }
    __host__ __device__ inline float g_theta_theta(float r) const { return r * r; }
    __host__ __device__ inline float g_phi_phi(float r, float theta) const {
        float s = sinf(theta); return r * r * s * s;
    }
};

// 4D position and momentum for geodesic integration
struct GeodesicState {
    float t, r, theta, phi;        // coordinates
    float pt, pr, ptheta, pphi;    // covariant momenta

    __host__ __device__ GeodesicState()
        : t(0), r(0), theta(0), phi(0), pt(0), pr(0), ptheta(0), pphi(0) {}
};

// Host-side lookup table (not used on device)
struct GeodesicLookupTable {
    float* impact_parameters;    // size num_entries
    float* deflection_angles;    // size num_entries
    float* bending_angles;       // size num_entries (NaN for non-escaped rays)
    float* closest_approaches;   // size num_entries
    int*   ray_status;           // size num_entries (compat with existing integer codes)
    int    num_entries;
    float  mass;
};

#endif // BH_TYPES_H
