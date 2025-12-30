/**
 * @file bh_schwarzschild.h
 * @brief Schwarzschild geodesic integration using the Hamiltonian formulation.
 *
 * This header provides functions for integrating null geodesics (photon paths)
 * in Schwarzschild spacetime using the Hamiltonian formulation with a 4th-order
 * Runge-Kutta (RK4) integrator.
 *
 * The Hamiltonian for a free particle in curved spacetime is:
 *   H = (1/2) g^{mu nu} p_mu p_nu
 *
 * For null geodesics (photons), H = 0. Hamilton's equations give:
 *   dx^mu/d(lambda) = dH/dp_mu = g^{mu nu} p_nu
 *   dp_mu/d(lambda) = -dH/dx^mu = -(1/2) (d g^{ab}/dx^mu) p_a p_b
 *
 * Key features:
 * - Hamiltonian constraint monitoring for accuracy validation
 * - Automatic singularity avoidance near horizon and poles
 * - Fixed-step RK4 integration with external step size control
 * - Device-only integration functions for CUDA kernel use
 *
 * @note Coordinates: Schwarzschild (t, r, theta, phi).
 * @note Metric signature: (-,+,+,+).
 * @note Units: Geometric units (G = c = 1).
 *
 * @see bh_common.h for constants and utilities
 * @see bh_types.h for GeodesicState and SchwarzschildMetric
 * @see bh_kerr.h for Kerr spacetime integration
 */

#ifndef BH_SCHWARZSCHILD_H
#define BH_SCHWARZSCHILD_H

#include "bh_common.h"
#include "bh_types.h"

/* ============================================================================
 * HAMILTONIAN CONSTRAINT
 * ============================================================================
 */

/**
 * @brief Compute the Hamiltonian for a geodesic state in Schwarzschild spacetime.
 *
 * Evaluates H = (1/2) g^{mu nu} p_mu p_nu using the inverse Schwarzschild metric:
 *   g^{tt} = -1/(1 - Rs/r)
 *   g^{rr} = (1 - Rs/r)
 *   g^{theta theta} = 1/r^2
 *   g^{phi phi} = 1/(r^2 sin^2(theta))
 *
 * For a null geodesic, H should be zero (or very small due to numerical error).
 * Monitoring |H| during integration provides a measure of accuracy.
 *
 * @param s The geodesic state containing position and covariant momenta.
 * @param Rs Schwarzschild radius Rs = 2M.
 * @return The Hamiltonian value H. Should be ~0 for null geodesics.
 *
 * @note Protected against division by zero at r = Rs and theta = 0, pi.
 *
 * Example usage:
 * @code
 *   float H = hamiltonian_schwarzschild(state, metric.Rs);
 *   if (fabsf(H) > 1e-6f) {
 *       // Significant Hamiltonian drift - accuracy degraded
 *   }
 * @endcode
 */
__host__ __device__ inline float hamiltonian_schwarzschild(const GeodesicState& s, float Rs) {
    float r = s.r;
    float th = s.theta;

    // Metric factor A = 1 - Rs/r, protected against horizon
    float A = fmaxf(BH_EPS, 1.0f - Rs / r);

    // Inverse metric components
    float gtt = -1.0f / A;
    float grr = A;
    float gthth = 1.0f / (r * r);
    float st2 = fmaxf(BH_EPS_SIN, sinf(th) * sinf(th));
    float gphph = 1.0f / (r * r * st2);

    // H = (1/2) g^{mu nu} p_mu p_nu
    float H = 0.5f * (gtt * s.pt * s.pt +
                      grr * s.pr * s.pr +
                      gthth * s.ptheta * s.ptheta +
                      gphph * s.pphi * s.pphi);
    return H;
}

/* ============================================================================
 * GEODESIC DERIVATIVES
 * ============================================================================
 */

/**
 * @brief Compute the derivatives of a geodesic state using Hamilton's equations.
 *
 * Given a geodesic state (position and covariant momenta), computes the
 * derivatives with respect to affine parameter lambda:
 *
 * Position derivatives (velocity):
 *   dt/d(lambda) = g^{tt} p_t
 *   dr/d(lambda) = g^{rr} p_r
 *   d(theta)/d(lambda) = g^{theta theta} p_theta
 *   d(phi)/d(lambda) = g^{phi phi} p_phi
 *
 * Momentum derivatives:
 *   dp_t/d(lambda) = 0  (t is cyclic - energy conserved)
 *   dp_r/d(lambda) = -(1/2) (d g^{ab}/dr) p_a p_b
 *   dp_theta/d(lambda) = -(1/2) (d g^{ab}/d(theta)) p_a p_b
 *   dp_phi/d(lambda) = 0  (phi is cyclic - angular momentum conserved)
 *
 * @param state Input geodesic state (position and momenta).
 * @param metric Schwarzschild metric (provides Rs).
 * @param[out] derivatives Output derivatives for RK4 integration.
 *
 * @note Device-only function (__device__).
 * @note If r <= Rs + 0.001 (near/inside horizon), returns zero derivatives
 *       to halt integration gracefully.
 *
 * @warning The derivatives of g^{phi phi} w.r.t. theta involve 1/sin^3(theta),
 *          which diverges at the poles. Protected with a floor value.
 */
__device__ inline void computeGeodesicDerivatives(const GeodesicState& state,
                                                  const SchwarzschildMetric& metric,
                                                  GeodesicState& derivatives) {
    float r = state.r;
    float theta = state.theta;
    float Rs = metric.Rs;

    // Near or inside event horizon - return zero derivatives to halt integration
    if (r <= Rs + 0.001f || r <= 0.001f) {
        derivatives = GeodesicState();
        return;
    }

    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float r2 = r * r;
    float r3 = r2 * r;
    float A = fmaxf(BH_EPS, 1.0f - Rs / r);  // g_rr^{-1} = g^{rr}, protected

    /* ------------------------------------------------------------------------
     * Position derivatives: dx^mu/d(lambda) = g^{mu nu} p_nu
     * ------------------------------------------------------------------------ */
    derivatives.t = (-1.0f / A) * state.pt;       // g^{tt} p_t
    derivatives.r = A * state.pr;                 // g^{rr} p_r
    derivatives.theta = state.ptheta / r2;        // g^{theta theta} p_theta

    // g^{phi phi} p_phi with pole protection
    float sin2 = fmaxf(BH_EPS_SIN, sin_theta * sin_theta);
    derivatives.phi = state.pphi / (r2 * sin2);

    /* ------------------------------------------------------------------------
     * Momentum derivatives: dp_mu/d(lambda) = -(1/2) (d g^{ab}/dx^mu) p_a p_b
     * ------------------------------------------------------------------------ */

    // dp_t/d(lambda) = 0 (t is cyclic - metric independent of t)
    derivatives.pt = 0.0f;

    // Radial derivatives of inverse metric components
    float dgtt_inv_dr   = Rs / (r2 * A * A);      // d(g^{tt})/dr
    float dgrr_inv_dr   = Rs / r2;                // d(g^{rr})/dr
    float dgthth_inv_dr = -2.0f / r3;             // d(g^{theta theta})/dr
    float dgphph_inv_dr = -2.0f / (r3 * sin2);    // d(g^{phi phi})/dr

    // dp_r/d(lambda) = -(1/2) sum over ab
    derivatives.pr = -0.5f * (
        dgtt_inv_dr   * state.pt     * state.pt +
        dgrr_inv_dr   * state.pr     * state.pr +
        dgthth_inv_dr * state.ptheta * state.ptheta +
        dgphph_inv_dr * state.pphi   * state.pphi
    );

    // dp_theta/d(lambda): only g^{phi phi} depends on theta
    // d(g^{phi phi})/d(theta) = d[1/(r^2 sin^2 theta)]/d(theta)
    //                        = -2 cos(theta) / (r^2 sin^3 theta)
    // dp_theta = -(1/2) * (-2 cos/sin^3) * r^{-2} * pphi^2
    //          = cos(theta) / (r^2 sin^3(theta)) * pphi^2
    float sin_theta_safe = (fabsf(sin_theta) < 1e-6f)
                           ? (sin_theta >= 0.0f ? 1e-6f : -1e-6f)
                           : sin_theta;
    derivatives.ptheta = (cos_theta / (r2 * sin2 * sin_theta_safe)) *
                         state.pphi * state.pphi;

    // dp_phi/d(lambda) = 0 (phi is cyclic - metric independent of phi)
    derivatives.pphi = 0.0f;
}

/* ============================================================================
 * RK4 INTEGRATOR
 * ============================================================================
 */

/**
 * @brief Perform one step of 4th-order Runge-Kutta integration for a geodesic.
 *
 * Advances the geodesic state by one step of size step_size in affine parameter.
 * Uses the classic RK4 method:
 *
 *   k1 = f(y_n)
 *   k2 = f(y_n + h/2 * k1)
 *   k3 = f(y_n + h/2 * k2)
 *   k4 = f(y_n + h * k3)
 *   y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
 *
 * The function also checks termination conditions and returns false if the
 * integration should stop.
 *
 * @param[in,out] state Geodesic state to advance. Modified in place.
 * @param metric Schwarzschild metric (provides Rs for horizon check).
 * @param step_size Step size in affine parameter lambda.
 * @return true if integration can continue, false if a termination condition
 *         is met (crossed horizon, hit pole, etc.).
 *
 * @note Device-only function (__device__).
 *
 * Termination conditions (returns false):
 * - r <= Rs + 0.001: Ray crossed or approached event horizon
 * - theta <= 0.001 or theta >= pi - 0.001: Ray hit coordinate pole
 *
 * Example usage:
 * @code
 *   SchwarzschildMetric metric(mass);
 *   GeodesicState state = initializeRay(...);
 *   float h = 0.01f;
 *   for (int step = 0; step < MAX_STEPS; ++step) {
 *       float scale = fminf(100.0f, fmaxf(1.0f, state.r / (0.5f * metric.Rs)));
 *       if (!integrateGeodesic(state, metric, h * scale)) {
 *           break;  // Ray terminated
 *       }
 *   }
 * @endcode
 */
__device__ inline bool integrateGeodesic(GeodesicState& state,
                                         const SchwarzschildMetric& metric,
                                         float step_size) {
    GeodesicState k1, k2, k3, k4;
    GeodesicState temp_state;

    // k1 = f(state)
    computeGeodesicDerivatives(state, metric, k1);

    // k2 = f(state + 0.5 * h * k1)
    temp_state.t = state.t + 0.5f * step_size * k1.t;
    temp_state.r = state.r + 0.5f * step_size * k1.r;
    temp_state.theta = state.theta + 0.5f * step_size * k1.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k1.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k1.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k1.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k1.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k1.pphi;

    computeGeodesicDerivatives(temp_state, metric, k2);

    // k3 = f(state + 0.5 * h * k2)
    temp_state.t = state.t + 0.5f * step_size * k2.t;
    temp_state.r = state.r + 0.5f * step_size * k2.r;
    temp_state.theta = state.theta + 0.5f * step_size * k2.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k2.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k2.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k2.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k2.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k2.pphi;

    computeGeodesicDerivatives(temp_state, metric, k3);

    // k4 = f(state + h * k3)
    temp_state.t = state.t + step_size * k3.t;
    temp_state.r = state.r + step_size * k3.r;
    temp_state.theta = state.theta + step_size * k3.theta;
    temp_state.phi = state.phi + step_size * k3.phi;
    temp_state.pt = state.pt + step_size * k3.pt;
    temp_state.pr = state.pr + step_size * k3.pr;
    temp_state.ptheta = state.ptheta + step_size * k3.ptheta;
    temp_state.pphi = state.pphi + step_size * k3.pphi;

    computeGeodesicDerivatives(temp_state, metric, k4);

    // Update state: y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    float h_6 = step_size / 6.0f;
    state.t      += h_6 * (k1.t      + 2.0f * k2.t      + 2.0f * k3.t      + k4.t);
    state.r      += h_6 * (k1.r      + 2.0f * k2.r      + 2.0f * k3.r      + k4.r);
    state.theta  += h_6 * (k1.theta  + 2.0f * k2.theta  + 2.0f * k3.theta  + k4.theta);
    state.phi    += h_6 * (k1.phi    + 2.0f * k2.phi    + 2.0f * k3.phi    + k4.phi);
    state.pt     += h_6 * (k1.pt     + 2.0f * k2.pt     + 2.0f * k3.pt     + k4.pt);
    state.pr     += h_6 * (k1.pr     + 2.0f * k2.pr     + 2.0f * k3.pr     + k4.pr);
    state.ptheta += h_6 * (k1.ptheta + 2.0f * k2.ptheta + 2.0f * k3.ptheta + k4.ptheta);
    state.pphi   += h_6 * (k1.pphi   + 2.0f * k2.pphi   + 2.0f * k3.pphi   + k4.pphi);

    // Check termination conditions
    // Returns true if integration can continue, false if it should stop
    return (state.r > metric.Rs + 0.001f &&
            state.theta > 0.001f && state.theta < BH_PI - 0.001f);
}

#endif // BH_SCHWARZSCHILD_H
