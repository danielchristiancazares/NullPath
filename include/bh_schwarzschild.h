// Header-only Schwarzschild geodesic math and RK4 stepper
#ifndef BH_SCHWARZSCHILD_H
#define BH_SCHWARZSCHILD_H

#include "bh_common.h"
#include "bh_types.h"

// Hamiltonian H = 1/2 g^{μν} p_μ p_ν (null geodesics => H≈0)
__host__ __device__ inline float hamiltonian_schwarzschild(const GeodesicState& s, float Rs) {
    float r = s.r;
    float th = s.theta;
    float A = fmaxf(BH_EPS, 1.0f - Rs / r);
    float gtt = -1.0f / A;
    float grr = A;
    float gthth = 1.0f / (r*r);
    float st2 = fmaxf(BH_EPS_SIN, sinf(th)*sinf(th));
    float gphph = 1.0f / (r*r*st2);
    float H = 0.5f * (gtt*s.pt*s.pt + grr*s.pr*s.pr + gthth*s.ptheta*s.ptheta + gphph*s.pphi*s.pphi);
    return H;
}

// Geodesic derivatives (Hamiltonian form) in Schwarzschild spacetime
__device__ inline void computeGeodesicDerivatives(const GeodesicState& state,
                                                  const SchwarzschildMetric& metric,
                                                  GeodesicState& derivatives) {
    float r = state.r;
    float theta = state.theta;
    float Rs = metric.Rs;

    // Avoid singularities
    if (r <= Rs + 0.001f || r <= 0.001f) {
        derivatives = GeodesicState();
        return;
    }

    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float r2 = r * r;
    float r3 = r2 * r;
    float A = 1.0f - Rs / r; // g_rr^{-1} = g^{rr}

    // dx^μ/dλ = g^{μν} p_ν
    derivatives.t = (-1.0f / A) * state.pt;       // g^{tt} p_t
    derivatives.r = A * state.pr;                 // g^{rr} p_r
    derivatives.theta = state.ptheta / r2;        // g^{θθ} p_θ
    derivatives.phi = state.pphi / (r2 * fmaxf(BH_EPS_SIN, sin_theta * sin_theta)); // g^{φφ} p_φ

    derivatives.pt = 0.0f; // t is cyclic

    // ∂_r g^{αβ}
    float sin2 = fmaxf(BH_EPS_SIN, sin_theta * sin_theta);
    float dgtt_inv_dr   = Rs / (r2 * A * A);
    float dgrr_inv_dr   = Rs / (r2);
    float dgthth_inv_dr = -2.0f / (r3);
    float dgphph_inv_dr = -2.0f / (r3 * sin2);

    derivatives.pr = -0.5f * (
        dgtt_inv_dr   * state.pt     * state.pt +
        dgrr_inv_dr   * state.pr     * state.pr +
        dgthth_inv_dr * state.ptheta * state.ptheta +
        dgphph_inv_dr * state.pphi   * state.pphi
    );

    // θ-momentum: only g^{φφ} depends on θ; protect sin^3 θ
    float sin_theta_safe = (fabsf(sin_theta) < 1e-6f) ? (sin_theta >= 0.0f ? 1e-6f : -1e-6f) : sin_theta;
    derivatives.ptheta = (cos_theta / (r2 * sin2 * sin_theta_safe)) * state.pphi * state.pphi;

    derivatives.pphi = 0.0f; // φ is cyclic
}

// 4th order Runge-Kutta integration for geodesics
__device__ inline bool integrateGeodesic(GeodesicState& state,
                                         const SchwarzschildMetric& metric,
                                         float step_size) {
    GeodesicState k1, k2, k3, k4;
    GeodesicState temp_state;

    computeGeodesicDerivatives(state, metric, k1);

    temp_state.t = state.t + 0.5f * step_size * k1.t;
    temp_state.r = state.r + 0.5f * step_size * k1.r;
    temp_state.theta = state.theta + 0.5f * step_size * k1.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k1.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k1.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k1.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k1.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k1.pphi;

    computeGeodesicDerivatives(temp_state, metric, k2);

    temp_state.t = state.t + 0.5f * step_size * k2.t;
    temp_state.r = state.r + 0.5f * step_size * k2.r;
    temp_state.theta = state.theta + 0.5f * step_size * k2.theta;
    temp_state.phi = state.phi + 0.5f * step_size * k2.phi;
    temp_state.pt = state.pt + 0.5f * step_size * k2.pt;
    temp_state.pr = state.pr + 0.5f * step_size * k2.pr;
    temp_state.ptheta = state.ptheta + 0.5f * step_size * k2.ptheta;
    temp_state.pphi = state.pphi + 0.5f * step_size * k2.pphi;

    computeGeodesicDerivatives(temp_state, metric, k3);

    temp_state.t = state.t + step_size * k3.t;
    temp_state.r = state.r + step_size * k3.r;
    temp_state.theta = state.theta + step_size * k3.theta;
    temp_state.phi = state.phi + step_size * k3.phi;
    temp_state.pt = state.pt + step_size * k3.pt;
    temp_state.pr = state.pr + step_size * k3.pr;
    temp_state.ptheta = state.ptheta + step_size * k3.ptheta;
    temp_state.pphi = state.pphi + step_size * k3.pphi;

    computeGeodesicDerivatives(temp_state, metric, k4);

    float h_6 = step_size / 6.0f;
    state.t      += h_6 * (k1.t      + 2.0f*k2.t      + 2.0f*k3.t      + k4.t);
    state.r      += h_6 * (k1.r      + 2.0f*k2.r      + 2.0f*k3.r      + k4.r);
    state.theta  += h_6 * (k1.theta  + 2.0f*k2.theta  + 2.0f*k3.theta  + k4.theta);
    state.phi    += h_6 * (k1.phi    + 2.0f*k2.phi    + 2.0f*k3.phi    + k4.phi);
    state.pt     += h_6 * (k1.pt     + 2.0f*k2.pt     + 2.0f*k3.pt     + k4.pt);
    state.pr     += h_6 * (k1.pr     + 2.0f*k2.pr     + 2.0f*k3.pr     + k4.pr);
    state.ptheta += h_6 * (k1.ptheta + 2.0f*k2.ptheta + 2.0f*k3.ptheta + k4.ptheta);
    state.pphi   += h_6 * (k1.pphi   + 2.0f*k2.pphi   + 2.0f*k3.pphi   + k4.pphi);

    // Termination bounds
    return (state.r > metric.Rs + 0.001f &&
            state.theta > 0.001f && state.theta < BH_PI - 0.001f);
}

#endif // BH_SCHWARZSCHILD_H

