// Minimal common header for shared constants and helpers.
// Non-invasive: nothing includes this yet; behavior unchanged.

#ifndef BH_COMMON_H
#define BH_COMMON_H

#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>

// Naming follows project conventions: constants in UPPER_SNAKE_CASE,
// types PascalCase, functions in snake_case.

// Math constants and epsilons
constexpr float  BH_PI   = 3.14159265358979323846f;
constexpr double BH_PI_D = 3.141592653589793238462643383279502884;
constexpr float BH_EPS = 1e-12f;
constexpr float BH_EPS_SIN = 1e-12f;

// Physical multipliers in geometric units (G=c=1)
constexpr float BH_SCHWARZSCHILD_MULTIPLIER = 2.0f;  // Rs = 2M
constexpr float BH_PHOTON_SPHERE_MULTIPLIER = 1.5f;  // r_ph = 1.5 Rs

// Integrator defaults (keep in sync with existing sources until migrated)
constexpr int   BH_MAX_INTEGRATION_STEPS_DEFAULT = 10000;
constexpr float BH_INTEGRATION_STEP_DEFAULT = 0.01f;
// Accurate-mode RK tolerances (used by FP64 adaptive integrator)
constexpr double BH_RK_REL_TOL   = 1e-10;
constexpr double BH_RK_ABS_TOL   = 1e-12;
constexpr double BH_RK_SAFETY    = 0.9;
constexpr double BH_RK_HMIN      = 1e-6;
constexpr double BH_RK_HMAX      = 0.2;

// Compatibility macros: keep existing code working without churn.
#ifndef SCHWARZSCHILD_MULTIPLIER
#define SCHWARZSCHILD_MULTIPLIER BH_SCHWARZSCHILD_MULTIPLIER
#endif
#ifndef PHOTON_SPHERE_MULTIPLIER
#define PHOTON_SPHERE_MULTIPLIER BH_PHOTON_SPHERE_MULTIPLIER
#endif
#ifndef MAX_INTEGRATION_STEPS
#define MAX_INTEGRATION_STEPS BH_MAX_INTEGRATION_STEPS_DEFAULT
#endif
#ifndef INTEGRATION_STEP_SIZE
#define INTEGRATION_STEP_SIZE BH_INTEGRATION_STEP_DEFAULT
#endif

// Provide PI_F for code that expects it
#ifndef PI_F
#define PI_F BH_PI
#endif

// Ray classification (mirrors current integer codes)
enum class RayStatus : int { Escaped = 0, Captured = 1, PhotonSphere = 2, Inconclusive = 3 };

// CUDA error checking helper (use for new code; existing code remains unchanged)
#define BH_CUDA_CHECK(expr)                                                     \
    do {                                                                        \
        cudaError_t _err = (expr);                                              \
        if (_err != cudaSuccess) {                                              \
            std::fprintf(stderr, "[CUDA] %s failed at %s:%d: %s\n",             \
                         #expr, __FILE__, __LINE__, cudaGetErrorString(_err));  \
            std::abort();                                                       \
        }                                                                       \
    } while (0)

// Small math helpers
__host__ __device__ inline float bh_clamp(float x, float lo, float hi) {
    return fminf(fmaxf(x, lo), hi);
}

__host__ __device__ inline float bh_clamp01(float x) {
    return bh_clamp(x, 0.0f, 1.0f);
}

__host__ __device__ inline float bh_safe_sin2(float theta) {
    float s = sinf(theta);
    return fmaxf(s * s, BH_EPS_SIN);
}

__host__ __device__ inline float bh_safe_div(float num, float den, float eps = BH_EPS) {
    float d = (fabsf(den) < eps) ? (den >= 0.0f ? eps : -eps) : den;
    return num / d;
}

#endif // BH_COMMON_H
