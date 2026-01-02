/**
 * @file bh_common.h
 * @brief Common constants, utilities, and macros for the NullPath geodesic tracer.
 *
 * This header provides the foundational layer for all NullPath modules:
 * - Physical and mathematical constants in geometric units (G=c=1)
 * - Numerical tolerances and integrator parameters
 * - Utility functions for safe arithmetic and clamping
 * - Ray classification enumeration
 * - CUDA error handling macros
 *
 * All other NullPath headers depend on this file. It is designed to be
 * header-only with inline functions suitable for both host and device code.
 *
 * @note Units: Geometric units where G = c = 1. Mass M sets the length scale.
 *       Schwarzschild radius Rs = 2M, photon sphere r_ph = 3M = 1.5*Rs.
 *
 * @see bh_types.h for data structures
 * @see bh_schwarzschild.h for Schwarzschild geodesics
 * @see bh_kerr.h for Kerr geodesics
 */

#ifndef BH_COMMON_H
#define BH_COMMON_H

#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>

/* ============================================================================
 * NAMING CONVENTIONS
 * ============================================================================
 * - Constants: UPPER_SNAKE_CASE with BH_ prefix
 * - Types/Structs: PascalCase
 * - Functions: snake_case with bh_ prefix for utilities
 * - Macros: UPPER_SNAKE_CASE with BH_ prefix
 * ============================================================================
 */

/* ============================================================================
 * MATHEMATICAL CONSTANTS
 * ============================================================================
 */

/**
 * @brief Pi in single precision.
 *
 * Used for angular calculations, coordinate transforms, and trigonometry.
 */
constexpr float BH_PI = 3.14159265358979323846f;

/**
 * @brief Pi in double precision.
 *
 * Used in FP64 geodesic integration paths for maximum accuracy.
 */
constexpr double BH_PI_D = 3.141592653589793238462643383279502884;

/**
 * @brief General-purpose small epsilon for avoiding division by zero.
 *
 * Used as a floor value when denominators approach zero.
 * Typical use: `denominator = fmax(actual_denom, BH_EPS)`
 *
 * @note Set to 1e-6f to be within practical FP32 precision (~7 decimal digits).
 *       Values smaller than ~1e-7 relative to the operands cannot be reliably
 *       distinguished in float arithmetic.
 */
constexpr float BH_EPS = 1e-6f;

/**
 * @brief Double-precision epsilon for FP64 code paths.
 *
 * Used in double-precision metric and integration functions.
 * Set to 1e-12 to be within practical FP64 precision (~15 decimal digits).
 */
constexpr double BH_EPS_D = 1e-12;

/**
 * @brief Epsilon for sin^2(theta) to avoid singularities at poles.
 *
 * The metric components g^{phi,phi} contain sin^2(theta) in the denominator,
 * which diverges at theta = 0 or pi. This epsilon prevents division by zero
 * while maintaining numerical accuracy away from the poles.
 *
 * @note Set to 1e-6f to be within practical FP32 precision.
 */
constexpr float BH_EPS_SIN = 1e-6f;

/**
 * @brief Double-precision epsilon for sin^2(theta) in FP64 code paths.
 */
constexpr double BH_EPS_SIN_D = 1e-12;

/* ============================================================================
 * PHYSICAL CONSTANTS (Geometric Units: G = c = 1)
 * ============================================================================
 */

/**
 * @brief Schwarzschild radius multiplier: Rs = 2M.
 *
 * In geometric units, the Schwarzschild radius (event horizon for a
 * non-rotating black hole) is exactly twice the mass parameter.
 */
constexpr float BH_SCHWARZSCHILD_MULTIPLIER = 2.0f;

/**
 * @brief Photon sphere radius as a multiple of Rs: r_ph = 1.5 * Rs = 3M.
 *
 * The photon sphere is the unstable circular orbit radius for photons.
 * Light rays with impact parameter b = b_crit asymptotically approach
 * this radius.
 */
constexpr float BH_PHOTON_SPHERE_MULTIPLIER = 1.5f;

/* ============================================================================
 * UNIT CONVERSION HELPERS
 * ============================================================================
 */

/**
 * @brief Compute Schwarzschild radius from mass.
 * @param M Black hole mass in geometric units.
 * @return Rs = 2M (event horizon radius).
 */
__host__ __device__ inline float bh_r_s(float M) {
    return BH_SCHWARZSCHILD_MULTIPLIER * M;
}

/**
 * @brief Compute photon sphere radius from Schwarzschild radius.
 * @param Rs Schwarzschild radius (Rs = 2M).
 * @return r_ph = 1.5 * Rs = 3M.
 */
__host__ __device__ inline float bh_r_ph(float Rs) {
    return BH_PHOTON_SPHERE_MULTIPLIER * Rs;
}

/**
 * @brief Compute critical impact parameter from mass.
 *
 * The critical impact parameter is the value of b = L/E at which a photon
 * asymptotically spirals into the photon sphere. For Schwarzschild:
 *   b_crit = 3*sqrt(3)*M
 *
 * @param M Black hole mass in geometric units.
 * @return b_crit = 3*sqrt(3)*M (approximately 5.196*M).
 */
__host__ __device__ inline float bh_bcrit_from_M(float M) {
    return 3.0f * sqrtf(3.0f) * M;
}

/**
 * @brief Compute critical impact parameter from Schwarzschild radius.
 *
 * Equivalent to bh_bcrit_from_M but takes Rs as input:
 *   b_crit = (3*sqrt(3)/2) * Rs
 *
 * @param Rs Schwarzschild radius (Rs = 2M).
 * @return b_crit = (3*sqrt(3)/2)*Rs (approximately 2.598*Rs).
 */
__host__ __device__ inline float bh_bcrit_from_Rs(float Rs) {
    return 0.5f * 3.0f * sqrtf(3.0f) * Rs;
}

/**
 * @brief Convert a length to units of M (mass).
 * @param x Length in geometric units.
 * @param M Black hole mass.
 * @return x/M (length in units of M).
 *
 * @note Protected against division by zero and overflow.
 */
__host__ __device__ inline float bh_to_M_units(float x, float M) {
    return x / fmaxf(M, BH_EPS);
}

/**
 * @brief Convert a length to units of Rs (Schwarzschild radius).
 * @param x Length in geometric units.
 * @param Rs Schwarzschild radius.
 * @return x/Rs (length in units of Rs).
 *
 * @note Protected against division by zero and overflow.
 */
__host__ __device__ inline float bh_to_Rs_units(float x, float Rs) {
    return x / fmaxf(Rs, BH_EPS);
}

/* ============================================================================
 * INTEGRATOR PARAMETERS
 * ============================================================================
 */

/**
 * @brief Default maximum number of integration steps per geodesic.
 *
 * Limits computational cost for rays that neither escape nor get captured
 * within a reasonable integration time. Rays exceeding this are marked
 * as RayStatus::Inconclusive.
 */
constexpr int BH_MAX_INTEGRATION_STEPS_DEFAULT = 10000;

/**
 * @brief Default base step size for RK4 integration.
 *
 * This is the step size in affine parameter lambda when r ~ Rs. The actual
 * step is scaled by radius: h = h_base * min(max(1, r/(0.5*Rs)), 100).
 *
 * @note Units: affine parameter (dimensionless in geometric units).
 */
constexpr float BH_INTEGRATION_STEP_DEFAULT = 0.01f;

/**
 * @brief Relative tolerance for adaptive RK4 (FP64 path).
 *
 * The step is accepted if the estimated error is less than:
 *   abs_tol + rel_tol * max(|y_full|, |y_half|)
 */
constexpr double BH_RK_REL_TOL = 1e-10;

/**
 * @brief Absolute tolerance for adaptive RK4 (FP64 path).
 *
 * Provides a floor for error tolerance when the solution approaches zero.
 */
constexpr double BH_RK_ABS_TOL = 1e-12;

/**
 * @brief Safety factor for adaptive step size adjustment.
 *
 * After estimating the optimal step size from error, multiply by this
 * factor to provide a margin of safety: h_new = safety * h_opt.
 */
constexpr double BH_RK_SAFETY = 0.9;

/**
 * @brief Minimum allowed step size for adaptive integration.
 *
 * Prevents the integrator from taking infinitesimally small steps near
 * singularities or turning points. If the required step is smaller,
 * accept the step anyway and continue.
 */
constexpr double BH_RK_HMIN = 1e-6;

/**
 * @brief Maximum allowed step size for adaptive integration.
 *
 * Prevents excessively large steps that could skip over important features
 * like disk crossings or photon sphere approaches.
 */
constexpr double BH_RK_HMAX = 0.2;

/* ============================================================================
 * BACKWARD COMPATIBILITY MACROS
 * ============================================================================
 * These macros maintain compatibility with existing code that uses the
 * older naming conventions. New code should use the BH_ prefixed constants.
 * ============================================================================
 */

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

#ifndef PI_F
#define PI_F BH_PI
#endif

/* ============================================================================
 * RAY CLASSIFICATION
 * ============================================================================
 */

/**
 * @brief Classification of ray termination conditions.
 *
 * After integrating a geodesic, each ray is classified based on its
 * final state. This enum uses explicit integer values for compatibility
 * with CUDA kernels that store results in int arrays.
 */
enum class RayStatus : int {
    /**
     * @brief Ray escaped to infinity.
     *
     * The ray reached a large radius (r > r_start) with outgoing momentum
     * (p_r > 0), indicating it will continue to infinity.
     */
    Escaped = 0,

    /**
     * @brief Ray was captured by the black hole.
     *
     * The ray crossed the event horizon (r <= Rs + epsilon) and will
     * reach the singularity.
     */
    Captured = 1,

    /**
     * @brief Ray approached the photon sphere.
     *
     * The ray's closest approach was within a small tolerance of
     * r_ph = 1.5*Rs, indicating it nearly achieved an unstable orbit.
     * These rays typically have b very close to b_crit.
     */
    PhotonSphere = 2,

    /**
     * @brief Ray status could not be determined.
     *
     * The integration reached MAX_INTEGRATION_STEPS without the ray
     * clearly escaping or being captured. This may indicate:
     * - Impact parameter very close to b_crit (long integration time)
     * - Numerical issues near the photon sphere
     * - Insufficient step budget
     */
    Inconclusive = 3
};

/* ============================================================================
 * CUDA ERROR HANDLING
 * ============================================================================
 */

/**
 * @brief CUDA error checking macro with file and line reporting.
 *
 * Wraps a CUDA API call and checks for errors. On failure, prints an
 * error message to stderr including the expression, file, line number,
 * and CUDA error string, then aborts the program.
 *
 * @param expr A CUDA API call that returns cudaError_t.
 *
 * Example usage:
 * @code
 *   BH_CUDA_CHECK(cudaMalloc(&ptr, size));
 *   BH_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
 * @endcode
 */
#define BH_CUDA_CHECK(expr)                                                     \
    do {                                                                        \
        cudaError_t _err = (expr);                                              \
        if (_err != cudaSuccess) {                                              \
            std::fprintf(stderr, "[CUDA] %s failed at %s:%d: %s\n",             \
                         #expr, __FILE__, __LINE__, cudaGetErrorString(_err));  \
            std::abort();                                                       \
        }                                                                       \
    } while (0)

/* ============================================================================
 * MATH UTILITY FUNCTIONS
 * ============================================================================
 */

/**
 * @brief Absolute value for float.
 */
__host__ __device__ inline float bh_abs(float x) {
    return fabsf(x);
}

/**
 * @brief Absolute value for double.
 */
__host__ __device__ inline double bh_abs(double x) {
    return fabs(x);
}

/**
 * @brief Clamp a value to a specified range.
 *
 * @param x Value to clamp.
 * @param lo Lower bound (inclusive).
 * @param hi Upper bound (inclusive).
 * @return Clamped value in [lo, hi].
 */
__host__ __device__ inline float bh_clamp(float x, float lo, float hi) {
    return fminf(fmaxf(x, lo), hi);
}

/**
 * @brief Clamp a value to [0, 1].
 *
 * Commonly used for color component clamping before output.
 *
 * @param x Value to clamp.
 * @return Clamped value in [0, 1].
 */
__host__ __device__ inline float bh_clamp01(float x) {
    return bh_clamp(x, 0.0f, 1.0f);
}

/**
 * @brief Compute sin^2(theta) with protection against zero at poles.
 *
 * The inverse metric component g^{phi,phi} contains 1/sin^2(theta),
 * which diverges at theta = 0 or pi. This function returns
 * max(sin^2(theta), BH_EPS_SIN) to prevent division by zero while
 * maintaining accuracy away from the poles.
 *
 * @param theta Polar angle in radians.
 * @return max(sin^2(theta), BH_EPS_SIN).
 */
__host__ __device__ inline float bh_safe_sin2(float theta) {
    float s = sinf(theta);
    return fmaxf(s * s, BH_EPS_SIN);
}

/**
 * @brief Safe division with protection against zero denominator.
 *
 * If the absolute value of the denominator is less than eps, replaces it
 * with +/-eps (preserving sign) before dividing.
 *
 * @param num Numerator.
 * @param den Denominator.
 * @param eps Minimum absolute value for denominator (default: BH_EPS).
 * @return num / safe_den where |safe_den| >= eps.
 *
 * @note Uses copysignf to correctly handle signed zero (-0.0f).
 */
__host__ __device__ inline float bh_safe_div(float num, float den, float eps = BH_EPS) {
    float d = (fabsf(den) < eps) ? copysignf(eps, den) : den;
    return num / d;
}

/* ============================================================================
 * COMPENSATED SUMMATION
 * ============================================================================
 */

/**
 * @brief Iterative Kahan-Babuska (Klein) compensated summation.
 *
 * Tracks first- and second-order compensation terms to reduce cancellation
 * in mixed-sign sums without resorting to double precision.
 */
template <typename Real>
struct BhKleinSum {
    Real sum;
    Real cs;
    Real ccs;

    __host__ __device__ inline BhKleinSum() : sum(Real(0)), cs(Real(0)), ccs(Real(0)) {}

    __host__ __device__ inline void add(Real x) {
        Real t = sum + x;
        Real c;
        if (bh_abs(sum) >= bh_abs(x)) {
            c = (sum - t) + x;
        } else {
            c = (x - t) + sum;
        }
        sum = t;

        t = cs + c;
        Real cc;
        if (bh_abs(cs) >= bh_abs(c)) {
            cc = (cs - t) + c;
        } else {
            cc = (c - t) + cs;
        }
        cs = t;
        ccs += cc;
    }

    __host__ __device__ inline Real result() const {
        return sum + cs + ccs;
    }
};

template <typename Real>
__host__ __device__ inline Real bh_rk4_sum(Real k1, Real k2, Real k3, Real k4) {
    BhKleinSum<Real> sum;
    sum.add(k1);
    sum.add(Real(2) * k2);
    sum.add(Real(2) * k3);
    sum.add(k4);
    return sum.result();
}

#endif // BH_COMMON_H
