/**
 * @file bh_kerr.h
 * @brief Kerr black hole geodesic integration, camera models, and Doppler calculations.
 *
 * This header provides comprehensive support for null geodesic integration in
 * Kerr spacetime (rotating black holes) using Boyer-Lindquist coordinates.
 *
 * Features:
 * - Covariant and contravariant Kerr metric components
 * - Closed-form analytic derivatives of inverse metric (no finite differences)
 * - Both FP32 and FP64 geodesic integration paths
 * - Fixed-step and adaptive RK4 integrators
 * - Static observer and ZAMO (Zero Angular Momentum Observer) camera models
 * - Doppler g-factor for disk emission with relativistic beaming
 * - Ergosphere detection for camera placement validation
 *
 * The Kerr metric in Boyer-Lindquist coordinates (t, r, theta, phi) is:
 *
 *   ds^2 = -(1 - 2Mr/Sigma) dt^2 - (4Mar sin^2(theta)/Sigma) dt dphi
 *          + (Sigma/Delta) dr^2 + Sigma dtheta^2
 *          + (A sin^2(theta)/Sigma) dphi^2
 *
 * where:
 *   Sigma = r^2 + a^2 cos^2(theta)
 *   Delta = r^2 - 2Mr + a^2
 *   A = (r^2 + a^2)^2 - a^2 Delta sin^2(theta)
 *   a = J/M (spin parameter, |a| <= M for physical black holes)
 *
 * Key radii:
 *   r_+ = M + sqrt(M^2 - a^2)  (outer event horizon)
 *   r_- = M - sqrt(M^2 - a^2)  (inner Cauchy horizon)
 *   r_ergo = M + sqrt(M^2 - a^2 cos^2(theta))  (ergosphere outer boundary)
 *
 * @note Coordinates: Boyer-Lindquist (t, r, theta, phi).
 * @note Metric signature: (-,+,+,+).
 * @note Units: Geometric units (G = c = 1).
 * @note Spin convention: a > 0 for prograde rotation (black hole spins counterclockwise
 *       when viewed from above the north pole).
 *
 * @see bh_common.h for constants and utilities
 * @see bh_types.h for GeodesicState
 * @see bh_schwarzschild.h for the a=0 limit
 */

#ifndef BH_KERR_H
#define BH_KERR_H

#include "bh_common.h"
#include "bh_types.h"

/* ============================================================================
 * KERR METRIC COMPONENTS (FP32)
 * ============================================================================
 */

/**
 * @brief Compute covariant Kerr metric components in Boyer-Lindquist coordinates.
 *
 * Returns the non-zero covariant metric components g_{mu nu}. The Kerr metric
 * has an off-diagonal g_{t phi} term encoding frame dragging.
 *
 * Components returned:
 *   g_tt = -(1 - 2Mr/Sigma)
 *   g_{t phi} = -2Mar sin^2(theta)/Sigma
 *   g_{phi phi} = A sin^2(theta)/Sigma
 *   g_rr = Sigma/Delta
 *   g_{theta theta} = Sigma
 *
 * @param M Black hole mass.
 * @param a Spin parameter (a = J/M, |a| <= M).
 * @param r Radial coordinate.
 * @param th Polar angle theta in radians.
 * @param[out] g_tt Covariant tt component.
 * @param[out] g_tphi Covariant t-phi component (off-diagonal).
 * @param[out] g_phiphi Covariant phi-phi component.
 * @param[out] g_rr Covariant rr component.
 * @param[out] g_thth Covariant theta-theta component.
 *
 * @note Protected against Delta = 0 (horizon) with BH_EPS floor.
 */
__device__ __host__ inline void kerr_blocks_cov(float M, float a, float r, float th,
                                               float& g_tt, float& g_tphi, float& g_phiphi,
                                               float& g_rr, float& g_thth) {
    float ct = cosf(th);
    float st = sinf(th);
    float st2 = st * st;
    float r2 = r * r;
    float a2 = a * a;
    float Sigma = r2 + a2 * ct * ct;
    float Delta = r2 - 2.0f * M * r + a2;
    float A = (r2 + a2) * (r2 + a2) - a2 * Delta * st2;

    g_tt     = -(1.0f - 2.0f * M * r / Sigma);
    g_tphi   = -2.0f * M * a * r * st2 / Sigma;
    g_phiphi = A * st2 / Sigma;
    g_rr     = Sigma / fmaxf(BH_EPS, Delta);
    g_thth   = Sigma;
}

/**
 * @brief Compute contravariant (inverse) Kerr metric components.
 *
 * Returns the non-zero contravariant metric components g^{mu nu} used in
 * the Hamiltonian formulation for raising indices.
 *
 * Components returned:
 *   g^{tt} = -A / (Sigma * Delta)
 *   g^{t phi} = -2Mar / (Sigma * Delta)
 *   g^{phi phi} = (Delta - a^2 sin^2(theta)) / (Sigma * Delta * sin^2(theta))
 *   g^{rr} = Delta / Sigma
 *   g^{theta theta} = 1 / Sigma
 *
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param r Radial coordinate.
 * @param th Polar angle theta.
 * @param[out] gtt Contravariant tt component.
 * @param[out] gtphi Contravariant t-phi component.
 * @param[out] gphiphi Contravariant phi-phi component.
 * @param[out] grr Contravariant rr component.
 * @param[out] gthth Contravariant theta-theta component.
 *
 * @note Protected against singularities at horizon (Delta=0) and poles (sin(theta)=0).
 */
__device__ __host__ inline void kerr_blocks_contra(float M, float a, float r, float th,
                                                  float& gtt, float& gtphi, float& gphiphi,
                                                  float& grr, float& gthth) {
    float ct = cosf(th);
    float st = sinf(th);
    float st2 = st * st;
    float r2 = r * r;
    float a2 = a * a;
    float Sigma = r2 + a2 * ct * ct;
    float Delta = r2 - 2.0f * M * r + a2;
    float A = (r2 + a2) * (r2 + a2) - a2 * Delta * st2;

    grr     = fmaxf(BH_EPS, Delta) / fmaxf(BH_EPS, Sigma);
    gthth   = 1.0f / fmaxf(BH_EPS, Sigma);
    gtt     = -A / fmaxf(BH_EPS, Sigma * Delta);
    gtphi   = -2.0f * M * a * r / fmaxf(BH_EPS, Sigma * Delta);
    gphiphi = (Delta - a2 * st2) / fmaxf(BH_EPS, Sigma * Delta * st2 + BH_EPS);
}

/* ============================================================================
 * ANALYTIC METRIC DERIVATIVES (FP32)
 * ============================================================================
 */

/**
 * @brief Compute analytic partial derivatives of inverse metric w.r.t. r and theta.
 *
 * Provides closed-form expressions for d(g^{mu nu})/dr and d(g^{mu nu})/d(theta)
 * for all contravariant metric components. These are needed for the momentum
 * evolution equations in the Hamiltonian formulation.
 *
 * Using analytic derivatives (rather than finite differences) eliminates
 * numerical differentiation errors and allows larger integration steps.
 *
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param r Radial coordinate.
 * @param th Polar angle theta.
 * @param[out] gtt_r Derivative d(g^{tt})/dr.
 * @param[out] gtphi_r Derivative d(g^{t phi})/dr.
 * @param[out] gphiphi_r Derivative d(g^{phi phi})/dr.
 * @param[out] grr_r Derivative d(g^{rr})/dr.
 * @param[out] gthth_r Derivative d(g^{theta theta})/dr.
 * @param[out] gtt_th Derivative d(g^{tt})/d(theta).
 * @param[out] gtphi_th Derivative d(g^{t phi})/d(theta).
 * @param[out] gphiphi_th Derivative d(g^{phi phi})/d(theta).
 * @param[out] grr_th Derivative d(g^{rr})/d(theta).
 * @param[out] gthth_th Derivative d(g^{theta theta})/d(theta).
 *
 * @note Protected against division by zero with floor values.
 */
__device__ __host__ inline void kerr_blocks_contra_derivs(float M, float a, float r, float th,
                                                        float& gtt_r, float& gtphi_r, float& gphiphi_r,
                                                        float& grr_r, float& gthth_r,
                                                        float& gtt_th, float& gtphi_th, float& gphiphi_th,
                                                        float& grr_th, float& gthth_th) {
    float s = sinf(th), c = cosf(th);
    float s2 = fmaxf(BH_EPS_SIN, s * s);
    float r2 = r * r;
    float a2 = a * a;
    float B = r2 + a2;
    float Sigma = r2 + a2 * c * c;
    float Delta = r2 - 2.0f * M * r + a2;
    float A = (B * B) - a2 * Delta * s2;
    float D = Sigma * Delta;
    float D2 = fmaxf(1e-24f, D * D);
    float s2_theta = 2.0f * s * c;  // d(sin^2)/d(theta) = 2 sin cos

    // Partial derivatives of intermediate quantities
    float Sigma_r = 2.0f * r;
    float Sigma_th = -2.0f * a2 * s * c;
    float Delta_r = 2.0f * (r - M);
    float A_r = 4.0f * r * B - a2 * Delta_r * s2;
    float A_th = -a2 * Delta * s2_theta;
    float D_r = Sigma_r * Delta + Sigma * Delta_r;
    float D_th = Sigma_th * Delta;  // Delta independent of theta

    // g^{tt} = -A/D  =>  d/dx(-A/D) = -(A_x D - A D_x)/D^2
    gtt_r = -(A_r * D - A * D_r) / D2;
    gtt_th = -(A_th * D - A * D_th) / D2;

    // g^{t phi} = -2aMr/D
    float N_tphi = 2.0f * a * M * r;
    float N_tphi_r = 2.0f * a * M;
    gtphi_r = -(N_tphi_r * D - N_tphi * D_r) / D2;
    gtphi_th = -(0.0f * D - N_tphi * D_th) / D2;  // N_tphi independent of theta

    // g^{rr} = Delta/Sigma
    grr_r = (Delta_r * Sigma - Delta * Sigma_r) / (Sigma * Sigma);
    grr_th = -(Delta * Sigma_th) / (Sigma * Sigma);

    // g^{theta theta} = 1/Sigma
    gthth_r = -Sigma_r / (Sigma * Sigma);
    gthth_th = -Sigma_th / (Sigma * Sigma);

    // g^{phi phi} = (Delta - a^2 sin^2)/(Delta Sigma sin^2)
    float N = Delta - a2 * s2;
    float N_r = Delta_r;
    float N_th = -a2 * s2_theta;
    float Den = Delta * Sigma * s2;
    float Den_r = Delta_r * Sigma * s2 + Delta * Sigma_r * s2;
    float Den_th = Delta * Sigma_th * s2 + Delta * Sigma * s2_theta;
    float Den2 = fmaxf(1e-24f, Den * Den);
    gphiphi_r = (N_r * Den - N * Den_r) / Den2;
    gphiphi_th = (N_th * Den - N * Den_th) / Den2;
}

/* ============================================================================
 * HAMILTONIAN (FP32)
 * ============================================================================
 */

/**
 * @brief Compute the Hamiltonian for a geodesic state in Kerr spacetime.
 *
 * Evaluates H = (1/2) g^{mu nu} p_mu p_nu including the off-diagonal t-phi term:
 *
 *   H = (1/2) [g^{tt} p_t^2 + 2 g^{t phi} p_t p_phi + g^{phi phi} p_phi^2
 *              + g^{rr} p_r^2 + g^{theta theta} p_theta^2]
 *
 * For a null geodesic, H should be zero.
 *
 * @param s Geodesic state with position and covariant momenta.
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @return Hamiltonian value (should be ~0 for null geodesics).
 */
__host__ __device__ inline float hamiltonian_kerr(const GeodesicState& s, float M, float a) {
    float gtt, gtphi, gphiphi, grr, gthth;
    kerr_blocks_contra(M, a, s.r, s.theta, gtt, gtphi, gphiphi, grr, gthth);

    BhKleinSum<float> sum;
    sum.add(gtt * s.pt * s.pt);
    sum.add(2.0f * gtphi * s.pt * s.pphi);
    sum.add(gphiphi * s.pphi * s.pphi);
    sum.add(grr * s.pr * s.pr);
    sum.add(gthth * s.ptheta * s.ptheta);
    float H = 0.5f * sum.result();
    return H;
}

/* ============================================================================
 * GEODESIC DERIVATIVES AND INTEGRATOR (FP32)
 * ============================================================================
 */

/**
 * @brief Compute geodesic derivatives in Kerr spacetime using Hamilton's equations.
 *
 * Position derivatives (velocity):
 *   dt/d(lambda) = g^{tt} p_t + g^{t phi} p_phi
 *   dr/d(lambda) = g^{rr} p_r
 *   d(theta)/d(lambda) = g^{theta theta} p_theta
 *   d(phi)/d(lambda) = g^{t phi} p_t + g^{phi phi} p_phi
 *
 * Momentum derivatives:
 *   dp_t/d(lambda) = 0  (stationarity)
 *   dp_phi/d(lambda) = 0  (axisymmetry)
 *   dp_r/d(lambda) = -(1/2) sum_{ab} (d g^{ab}/dr) p_a p_b
 *   dp_theta/d(lambda) = -(1/2) sum_{ab} (d g^{ab}/d(theta)) p_a p_b
 *
 * @param st Input geodesic state.
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param[out] d Output derivatives.
 *
 * @note Device-only function.
 * @note Returns zero derivatives if r <= r_+ + 1e-3 (near horizon).
 */
__device__ inline void computeGeodesicDerivativesKerr(const GeodesicState& st, float M, float a,
                                                      GeodesicState& d) {
    float r = st.r;
    float th = st.theta;

    // Outer horizon radius r_+ = M + sqrt(M^2 - a^2)
    float disc = fmaxf(0.0f, M * M - a * a);
    float r_plus = M + sqrtf(disc);
    if (r <= r_plus + 1e-3f || r <= 1e-6f) {
        d = GeodesicState();
        return;
    }

    // Get metric and derivatives
    float gtt, gtphi, gphiphi, grr, gthth;
    kerr_blocks_contra(M, a, r, th, gtt, gtphi, gphiphi, grr, gthth);

    float gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r;
    float gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th;
    kerr_blocks_contra_derivs(M, a, r, th,
                              gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r,
                              gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th);

    // Position derivatives: dx^mu/d(lambda) = g^{mu nu} p_nu
    d.t     = gtt   * st.pt + gtphi   * st.pphi;
    d.r     = grr   * st.pr;
    d.theta = gthth * st.ptheta;
    d.phi   = gtphi * st.pt + gphiphi * st.pphi;

    // Momentum derivatives: dp_mu/d(lambda) = -(1/2) (d g^{ab}/dx^mu) p_a p_b
    d.pt = 0.0f;    // Stationarity (metric independent of t)
    d.pphi = 0.0f;  // Axisymmetry (metric independent of phi)

    BhKleinSum<float> sum_r;
    sum_r.add(gtt_r * st.pt * st.pt);
    sum_r.add(2.0f * gtphi_r * st.pt * st.pphi);
    sum_r.add(gphiphi_r * st.pphi * st.pphi);
    sum_r.add(grr_r * st.pr * st.pr);
    sum_r.add(gthth_r * st.ptheta * st.ptheta);
    float S_r = sum_r.result();

    BhKleinSum<float> sum_th;
    sum_th.add(gtt_th * st.pt * st.pt);
    sum_th.add(2.0f * gtphi_th * st.pt * st.pphi);
    sum_th.add(gphiphi_th * st.pphi * st.pphi);
    sum_th.add(grr_th * st.pr * st.pr);
    sum_th.add(gthth_th * st.ptheta * st.ptheta);
    float S_th = sum_th.result();

    d.pr = -0.5f * S_r;
    d.ptheta = -0.5f * S_th;
}

/**
 * @brief Perform one RK4 step for Kerr geodesic integration (FP32).
 *
 * @param[in,out] st Geodesic state to advance.
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param h Step size in affine parameter.
 * @return true if integration can continue, false if terminated (horizon/pole).
 *
 * @note Device-only function.
 */
__device__ inline bool integrateGeodesicKerr(GeodesicState& st, float M, float a, float h) {
    GeodesicState k1, k2, k3, k4, tmp;

    computeGeodesicDerivativesKerr(st, M, a, k1);

    tmp.t = st.t + 0.5f*h*k1.t; tmp.r = st.r + 0.5f*h*k1.r;
    tmp.theta = st.theta + 0.5f*h*k1.theta; tmp.phi = st.phi + 0.5f*h*k1.phi;
    tmp.pt = st.pt + 0.5f*h*k1.pt; tmp.pr = st.pr + 0.5f*h*k1.pr;
    tmp.ptheta = st.ptheta + 0.5f*h*k1.ptheta; tmp.pphi = st.pphi + 0.5f*h*k1.pphi;

    computeGeodesicDerivativesKerr(tmp, M, a, k2);

    tmp.t = st.t + 0.5f*h*k2.t; tmp.r = st.r + 0.5f*h*k2.r;
    tmp.theta = st.theta + 0.5f*h*k2.theta; tmp.phi = st.phi + 0.5f*h*k2.phi;
    tmp.pt = st.pt + 0.5f*h*k2.pt; tmp.pr = st.pr + 0.5f*h*k2.pr;
    tmp.ptheta = st.ptheta + 0.5f*h*k2.ptheta; tmp.pphi = st.pphi + 0.5f*h*k2.pphi;

    computeGeodesicDerivativesKerr(tmp, M, a, k3);

    tmp.t = st.t + h*k3.t; tmp.r = st.r + h*k3.r;
    tmp.theta = st.theta + h*k3.theta; tmp.phi = st.phi + h*k3.phi;
    tmp.pt = st.pt + h*k3.pt; tmp.pr = st.pr + h*k3.pr;
    tmp.ptheta = st.ptheta + h*k3.ptheta; tmp.pphi = st.pphi + h*k3.pphi;

    computeGeodesicDerivativesKerr(tmp, M, a, k4);

    float h6 = h / 6.0f;
    st.t += h6 * bh_rk4_sum(k1.t, k2.t, k3.t, k4.t);
    st.r += h6 * bh_rk4_sum(k1.r, k2.r, k3.r, k4.r);
    st.theta += h6 * bh_rk4_sum(k1.theta, k2.theta, k3.theta, k4.theta);
    st.phi += h6 * bh_rk4_sum(k1.phi, k2.phi, k3.phi, k4.phi);
    st.pt += h6 * bh_rk4_sum(k1.pt, k2.pt, k3.pt, k4.pt);
    st.pr += h6 * bh_rk4_sum(k1.pr, k2.pr, k3.pr, k4.pr);
    st.ptheta += h6 * bh_rk4_sum(k1.ptheta, k2.ptheta, k3.ptheta, k4.ptheta);
    st.pphi += h6 * bh_rk4_sum(k1.pphi, k2.pphi, k3.pphi, k4.pphi);

    // Termination check
    float disc = fmaxf(0.0f, M*M - a*a);
    float r_plus = M + sqrtf(disc);
    return (st.r > r_plus + 1e-3f && st.theta > 1e-3f && st.theta < (float)BH_PI - 1e-3f);
}

/* ============================================================================
 * CAMERA INITIALIZATION (FP32)
 * ============================================================================
 */

/**
 * @brief Initialize a ray from a static observer camera in Kerr spacetime.
 *
 * Constructs the initial geodesic state for a photon emitted from a camera
 * at position (r_cam, theta_cam, phi_cam) with local direction (n_r, n_theta, n_phi).
 *
 * The static observer has 4-velocity u^mu proportional to the timelike Killing
 * vector (1, 0, 0, 0). This observer is only valid outside the ergosphere
 * where g_{tt} < 0.
 *
 * An orthonormal tetrad is constructed and the local direction components
 * are transformed to covariant 4-momentum.
 *
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param r_cam Camera radial position.
 * @param theta_cam Camera polar angle.
 * @param phi_cam Camera azimuthal angle.
 * @param n_r Local radial direction component (typically -1 for ingoing).
 * @param n_theta Local polar direction component.
 * @param n_phi Local azimuthal direction component.
 * @return Initialized GeodesicState ready for integration.
 *
 * @warning Only valid outside the ergosphere. Use kerr_static_timelike() to check.
 * @see init_ray_from_camera_zamo_d for ZAMO camera (valid inside ergosphere).
 */
__device__ __host__ inline GeodesicState init_ray_from_camera_general(float M, float a,
                                                                     float r_cam,
                                                                     float theta_cam,
                                                                     float phi_cam,
                                                                     float n_r, float n_theta, float n_phi) {
    GeodesicState s;
    s.t = 0.0f; s.r = r_cam; s.theta = theta_cam; s.phi = phi_cam;

    float g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov(M, a, s.r, s.theta, g_tt, g_tphi, g_phiphi, g_rr, g_thth);

    // Static observer tetrad (only valid when g_tt < 0)
    float ut = rsqrtf(fmaxf(BH_EPS, -g_tt));
    float er_r = rsqrtf(fmaxf(BH_EPS, g_rr));
    float eth_th = rsqrtf(fmaxf(BH_EPS, g_thth));

    float ephi_t, ephi_phi;
    if (fabsf(g_tphi) > 1e-12f) {
        float beta = -g_tt / g_tphi;
        BhKleinSum<float> norm2_sum;
        norm2_sum.add(g_tt);
        norm2_sum.add(2.0f * g_tphi * beta);
        norm2_sum.add(g_phiphi * beta * beta);
        float norm2 = norm2_sum.result();
        float inv = rsqrtf(fabsf(norm2) + 1e-12f);
        ephi_t = inv;
        ephi_phi = beta * inv;
    } else {
        ephi_t = 0.0f;
        ephi_phi = rsqrtf(fmaxf(BH_EPS, g_phiphi));
    }

    // Local energy (normalized to 1)
    const float E_loc = 1.0f;

    // Contravariant 4-momentum in coordinate basis
    float pt_contra = E_loc * (ut + n_phi * ephi_t);
    float pr_contra = E_loc * (n_r * er_r);
    float pth_contra = E_loc * (n_theta * eth_th);
    float pph_contra = E_loc * (n_phi * ephi_phi);

    // Lower indices: p_mu = g_{mu nu} p^nu
    s.pt = g_tt * pt_contra + g_tphi * pph_contra;
    s.pr = g_rr * pr_contra;
    s.ptheta = g_thth * pth_contra;
    s.pphi = g_tphi * pt_contra + g_phiphi * pph_contra;

    return s;
}

/* ============================================================================
 * DOPPLER FACTOR (FP32)
 * ============================================================================
 */

/**
 * @brief Compute the Doppler g-factor for emission from a Keplerian disk.
 *
 * The g-factor relates observed and emitted specific intensities:
 *   I_obs = g^3 * I_em  (for optically thick emission)
 *   I_obs = g^4 * I_em  (for optically thin emission)
 *
 * This function computes g = E_obs / E_em where:
 * - E_obs = -p_t (photon energy at infinity)
 * - E_em = -k_mu u^mu (photon energy in emitter frame)
 * - u^mu is the 4-velocity of matter on a prograde circular Keplerian orbit
 *
 * @param st Geodesic state at the emission point (must be at theta = pi/2).
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @return Doppler factor g >= 0.
 *
 * @note Assumes equatorial emission (theta = pi/2).
 * @note Uses prograde Keplerian angular velocity Omega = 1/(a + r^{3/2}/M).
 */
__device__ inline float doppler_g_general(const GeodesicState& st, float M, float a) {
    float r = st.r;
    float th = (float)BH_PI * 0.5f;  // Evaluate at disk plane

    float g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov(M, a, r, th, g_tt, g_tphi, g_phiphi, g_rr, g_thth);

    // Prograde Keplerian angular velocity in BL coordinates
    float r_over_M = r / fmaxf(BH_EPS, M);
    float Omega = 1.0f / (a / fmaxf(BH_EPS, M) + powf(fmaxf(1.0f, r_over_M), 1.5f)) / fmaxf(BH_EPS, M);

    // Normalization: u_mu u^mu = -1 => u^t = 1/sqrt(-g_tt - 2*g_tphi*Omega - g_phiphi*Omega^2)
    BhKleinSum<float> denom_sum;
    denom_sum.add(g_tt);
    denom_sum.add(2.0f * g_tphi * Omega);
    denom_sum.add(g_phiphi * Omega * Omega);
    float denom = -denom_sum.result();
    float u_t = rsqrtf(fmaxf(BH_EPS, denom));
    float u_phi = Omega * u_t;

    // g = E_inf / E_em = -p_t / (-k_mu u^mu)
    float E_inf = -st.pt;
    BhKleinSum<float> k_u_sum;
    k_u_sum.add(st.pt * u_t);
    k_u_sum.add(st.pphi * u_phi);
    float k_u = -k_u_sum.result();
    float g = E_inf / fmaxf(1e-9f, k_u);

    return fmaxf(0.0f, g);
}

/* ============================================================================
 * DOUBLE PRECISION (FP64) SECTION
 * ============================================================================
 * The following provides double-precision versions of the geodesic integration
 * routines for high-accuracy work, validation, and rays near critical parameters.
 * ============================================================================
 */

#pragma region FP64_geodesics

/**
 * @brief Double-precision geodesic state.
 *
 * Equivalent to GeodesicState but using double precision for improved
 * accuracy in long integrations or near critical impact parameters.
 */
struct GeodesicStateD {
    double t, r, theta, phi;           ///< Position coordinates
    double pt, pr, ptheta, pphi;       ///< Covariant 4-momentum
};

/**
 * @brief Compute covariant Kerr metric components (double precision).
 *
 * Returns the non-zero covariant metric components g_{mu nu} in double precision.
 *
 * @param M Black hole mass.
 * @param a Spin parameter (a = J/M, |a| <= M).
 * @param r Radial coordinate.
 * @param th Polar angle theta in radians.
 * @param[out] g_tt Covariant tt component.
 * @param[out] g_tphi Covariant t-phi component (off-diagonal).
 * @param[out] g_phiphi Covariant phi-phi component.
 * @param[out] g_rr Covariant rr component.
 * @param[out] g_thth Covariant theta-theta component.
 *
 * @note Protected against Delta = 0 (horizon) with BH_EPS_D floor.
 */
__device__ __host__ inline void kerr_blocks_cov_d(double M, double a, double r, double th,
                                                  double& g_tt, double& g_tphi, double& g_phiphi,
                                                  double& g_rr, double& g_thth) {
    double ct = cos(th);
    double st = sin(th);
    double st2 = st * st;
    double r2 = r * r;
    double a2 = a * a;
    double Sigma = r2 + a2 * ct * ct;
    double Delta = r2 - 2.0 * M * r + a2;
    double A = (r2 + a2) * (r2 + a2) - a2 * Delta * st2;

    g_tt     = -(1.0 - 2.0 * M * r / fmax(BH_EPS_D, Sigma));
    g_tphi   = -2.0 * M * a * r * st2 / fmax(BH_EPS_D, Sigma);
    g_phiphi = A * st2 / fmax(BH_EPS_D, Sigma);
    g_rr     = Sigma / fmax(BH_EPS_D, Delta);
    g_thth   = Sigma;
}

/**
 * @brief Compute contravariant Kerr metric components (double precision).
 *
 * @see kerr_blocks_contra for parameter documentation.
 * @note Protected against singularities at horizon (Delta=0) and poles (sin(theta)=0).
 */
__device__ __host__ inline void kerr_blocks_contra_d(double M, double a, double r, double th,
                                                    double& gtt, double& gtphi, double& gphiphi,
                                                    double& grr, double& gthth) {
    double ct = cos(th);
    double st = sin(th);
    double st2 = fmax(BH_EPS_SIN_D, st * st);
    double r2 = r*r;
    double a2 = a*a;
    double Sigma = fmax(BH_EPS_D, r2 + a2 * ct*ct);
    double Delta = r2 - 2.0*M*r + a2;
    double A = (r2 + a2)*(r2 + a2) - a2 * Delta * st2;

    grr     = fmax(BH_EPS_D, Delta) / Sigma;
    gthth   = 1.0 / Sigma;
    gtt     = -A / fmax(BH_EPS_D, Sigma * Delta);
    gtphi   = -2.0 * M * a * r / fmax(BH_EPS_D, Sigma * Delta);
    gphiphi = (Delta - a2 * st2) / fmax(BH_EPS_D, Sigma * Delta * st2);
}

/**
 * @brief Compute analytic metric derivatives (double precision).
 *
 * @see kerr_blocks_contra_derivs for parameter documentation.
 */
__device__ __host__ inline void kerr_blocks_contra_derivs_d(double M, double a, double r, double th,
                                                           double& gtt_r, double& gtphi_r, double& gphiphi_r,
                                                           double& grr_r, double& gthth_r,
                                                           double& gtt_th, double& gtphi_th, double& gphiphi_th,
                                                           double& grr_th, double& gthth_th) {
    double s = sin(th), c = cos(th);
    double s2 = (s*s > 1e-30 ? s*s : 1e-30);
    double r2 = r*r;
    double a2 = a*a;
    double B = r2 + a2;
    double Sigma = r2 + a2 * c*c;
    double Delta = r2 - 2.0*M*r + a2;
    double A = (B*B) - a2 * Delta * s2;
    double D = Sigma * Delta;
    double D2 = D*D + 1e-48;
    double s2_th = 2.0 * s * c;

    double Sigma_r = 2.0 * r;
    double Sigma_th = -2.0 * a2 * s * c;
    double Delta_r = 2.0 * (r - M);
    double A_r = 4.0 * r * B - a2 * Delta_r * s2;
    double A_th = -a2 * Delta * s2_th;
    double D_r = Sigma_r * Delta + Sigma * Delta_r;
    double D_th = Sigma_th * Delta;

    gtt_r = -(A_r * D - A * D_r) / D2;
    gtt_th = -(A_th * D - A * D_th) / D2;

    double Nt = 2.0 * M * a * r;
    double Nt_r = 2.0 * M * a;
    gtphi_r = -(Nt_r * D - Nt * D_r) / D2;
    gtphi_th = -(0.0 * D - Nt * D_th) / D2;

    grr_r = (Delta_r * Sigma - Delta * Sigma_r) / (Sigma * Sigma);
    grr_th = -(Delta * Sigma_th) / (Sigma * Sigma);
    gthth_r = -Sigma_r / (Sigma * Sigma);
    gthth_th = -Sigma_th / (Sigma * Sigma);

    double N = Delta - a2 * s2;
    double N_r = Delta_r;
    double N_th = -a2 * s2_th;
    double Den = Delta * Sigma * s2 + 1e-48;
    double Den_r = Delta_r * Sigma * s2 + Delta * Sigma_r * s2;
    double Den_th = Delta * Sigma_th * s2 + Delta * Sigma * s2_th;
    double Den2 = Den * Den + 1e-48;

    gphiphi_r = (N_r * Den - N * Den_r) / Den2;
    gphiphi_th = (N_th * Den - N * Den_th) / Den2;
}

/**
 * @brief Compute geodesic derivatives (double precision).
 *
 * @see computeGeodesicDerivativesKerr for parameter documentation.
 */
__device__ inline void computeGeodesicDerivativesKerrD(const GeodesicStateD& st, double M, double a,
                                                       GeodesicStateD& d) {
    double r = st.r;
    double th = st.theta;

    double disc = M*M - a*a;
    if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    if (r <= r_plus + 1e-12 || r <= 1e-12) {
        d = GeodesicStateD();
        return;
    }

    double gtt, gtphi, gphiphi, grr, gthth;
    kerr_blocks_contra_d(M, a, r, th, gtt, gtphi, gphiphi, grr, gthth);

    double gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r;
    double gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th;
    kerr_blocks_contra_derivs_d(M, a, r, th,
                                gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r,
                                gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th);

    d.t     = gtt   * st.pt + gtphi   * st.pphi;
    d.r     = grr   * st.pr;
    d.theta = gthth * st.ptheta;
    d.phi   = gtphi * st.pt + gphiphi * st.pphi;

    d.pt = 0.0;
    d.pphi = 0.0;

    BhKleinSum<double> sum_r;
    sum_r.add(gtt_r * st.pt * st.pt);
    sum_r.add(2.0 * gtphi_r * st.pt * st.pphi);
    sum_r.add(gphiphi_r * st.pphi * st.pphi);
    sum_r.add(grr_r * st.pr * st.pr);
    sum_r.add(gthth_r * st.ptheta * st.ptheta);
    double S_r = sum_r.result();

    BhKleinSum<double> sum_th;
    sum_th.add(gtt_th * st.pt * st.pt);
    sum_th.add(2.0 * gtphi_th * st.pt * st.pphi);
    sum_th.add(gphiphi_th * st.pphi * st.pphi);
    sum_th.add(grr_th * st.pr * st.pr);
    sum_th.add(gthth_th * st.ptheta * st.ptheta);
    double S_th = sum_th.result();

    d.pr = -0.5 * S_r;
    d.ptheta = -0.5 * S_th;
}

/**
 * @brief Fixed-step RK4 integrator (double precision).
 *
 * @see integrateGeodesicKerr for parameter documentation.
 */
__device__ inline bool integrateGeodesicKerrD(GeodesicStateD& st, double M, double a, double h) {
    GeodesicStateD k1, k2, k3, k4, tmp;

    computeGeodesicDerivativesKerrD(st, M, a, k1);

    tmp.t = st.t + 0.5*h*k1.t; tmp.r = st.r + 0.5*h*k1.r;
    tmp.theta = st.theta + 0.5*h*k1.theta; tmp.phi = st.phi + 0.5*h*k1.phi;
    tmp.pt = st.pt + 0.5*h*k1.pt; tmp.pr = st.pr + 0.5*h*k1.pr;
    tmp.ptheta = st.ptheta + 0.5*h*k1.ptheta; tmp.pphi = st.pphi + 0.5*h*k1.pphi;

    computeGeodesicDerivativesKerrD(tmp, M, a, k2);

    tmp.t = st.t + 0.5*h*k2.t; tmp.r = st.r + 0.5*h*k2.r;
    tmp.theta = st.theta + 0.5*h*k2.theta; tmp.phi = st.phi + 0.5*h*k2.phi;
    tmp.pt = st.pt + 0.5*h*k2.pt; tmp.pr = st.pr + 0.5*h*k2.pr;
    tmp.ptheta = st.ptheta + 0.5*h*k2.ptheta; tmp.pphi = st.pphi + 0.5*h*k2.pphi;

    computeGeodesicDerivativesKerrD(tmp, M, a, k3);

    tmp.t = st.t + h*k3.t; tmp.r = st.r + h*k3.r;
    tmp.theta = st.theta + h*k3.theta; tmp.phi = st.phi + h*k3.phi;
    tmp.pt = st.pt + h*k3.pt; tmp.pr = st.pr + h*k3.pr;
    tmp.ptheta = st.ptheta + h*k3.ptheta; tmp.pphi = st.pphi + h*k3.pphi;

    computeGeodesicDerivativesKerrD(tmp, M, a, k4);

    double h6 = h / 6.0;
    st.t += h6 * bh_rk4_sum(k1.t, k2.t, k3.t, k4.t);
    st.r += h6 * bh_rk4_sum(k1.r, k2.r, k3.r, k4.r);
    st.theta += h6 * bh_rk4_sum(k1.theta, k2.theta, k3.theta, k4.theta);
    st.phi += h6 * bh_rk4_sum(k1.phi, k2.phi, k3.phi, k4.phi);
    st.pt += h6 * bh_rk4_sum(k1.pt, k2.pt, k3.pt, k4.pt);
    st.pr += h6 * bh_rk4_sum(k1.pr, k2.pr, k3.pr, k4.pr);
    st.ptheta += h6 * bh_rk4_sum(k1.ptheta, k2.ptheta, k3.ptheta, k4.ptheta);
    st.pphi += h6 * bh_rk4_sum(k1.pphi, k2.pphi, k3.pphi, k4.pphi);

    double disc = M*M - a*a;
    if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    return (st.r > r_plus + 1e-9 && st.theta > 1e-9 && st.theta < BH_PI_D - 1e-9);
}

/**
 * @brief Initialize ray from static observer camera (double precision).
 *
 * @see init_ray_from_camera_general for parameter documentation.
 */
__device__ __host__ inline GeodesicStateD init_ray_from_camera_general_d(double M, double a,
                                                                        double r_cam,
                                                                        double theta_cam,
                                                                        double phi_cam,
                                                                        double n_r, double n_theta, double n_phi) {
    GeodesicStateD s{};
    s.t = 0.0; s.r = r_cam; s.theta = theta_cam; s.phi = phi_cam;

    // Use true FP64 metric for tetrad construction
    double g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov_d(M, a, s.r, s.theta, g_tt, g_tphi, g_phiphi, g_rr, g_thth);

    double ut = rsqrt(fmax(BH_EPS_D, -g_tt));
    double er_r = rsqrt(fmax(BH_EPS_D, g_rr));
    double eth_th = rsqrt(fmax(BH_EPS_D, g_thth));

    // Match FP32 tetrad construction: special case when g_tphi is small (Schwarzschild limit)
    double ephi_t, ephi_phi;
    if (fabs(g_tphi) > BH_EPS_D) {
        double beta = -g_tt / g_tphi;
        BhKleinSum<double> norm2_sum;
        norm2_sum.add(g_tt);
        norm2_sum.add(2.0 * g_tphi * beta);
        norm2_sum.add(g_phiphi * beta * beta);
        double norm2 = norm2_sum.result();
        double inv = rsqrt(fmax(BH_EPS_D, fabs(norm2)));
        ephi_t = inv;
        ephi_phi = beta * inv;
    } else {
        ephi_t = 0.0;
        ephi_phi = rsqrt(fmax(BH_EPS_D, g_phiphi));
    }

    double E_loc = 1.0;
    double pt_contra = E_loc * (ut + n_phi * ephi_t);
    double pr_contra = E_loc * (n_r * er_r);
    double pth_contra = E_loc * (n_theta * eth_th);
    double pph_contra = E_loc * (n_phi * ephi_phi);

    s.pt = g_tt * pt_contra + g_tphi * pph_contra;
    s.pr = g_rr * pr_contra;
    s.ptheta = g_thth * pth_contra;
    s.pphi = g_tphi * pt_contra + g_phiphi * pph_contra;

    return s;
}

/**
 * @brief Compute Doppler g-factor (double precision).
 *
 * @see doppler_g_general for parameter documentation.
 */
__device__ inline double doppler_g_general_d(const GeodesicStateD& st, double M, double a) {
    double r = st.r;
    double th = BH_PI_D * 0.5;

    // Use true FP64 metric
    double g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov_d(M, a, r, th, g_tt, g_tphi, g_phiphi, g_rr, g_thth);

    // Prograde Keplerian angular velocity with consistent floor values
    double M_safe = fmax(BH_EPS_D, M);
    double r_over_M = r / M_safe;
    double Omega = 1.0 / (a / M_safe + pow(fmax(1.0, r_over_M), 1.5)) / M_safe;

    BhKleinSum<double> denom_sum;
    denom_sum.add(g_tt);
    denom_sum.add(2.0 * g_tphi * Omega);
    denom_sum.add(g_phiphi * Omega * Omega);
    double denom = -denom_sum.result();
    double u_t = rsqrt(denom > 1e-24 ? denom : 1e-24);
    double u_phi = Omega * u_t;

    double E_inf = -st.pt;
    BhKleinSum<double> k_u_sum;
    k_u_sum.add(st.pt * u_t);
    k_u_sum.add(st.pphi * u_phi);
    double k_u = -k_u_sum.result();
    double g = E_inf / (k_u > 1e-18 ? k_u : 1e-18);

    return g > 0.0 ? g : 0.0;
}

/**
 * @brief Check if a static observer is timelike at given position.
 *
 * A static observer (worldline along the timelike Killing vector) exists
 * only where g_{tt} < 0. Inside the ergosphere, g_{tt} >= 0 and static
 * observers cannot exist.
 *
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param r Radial coordinate.
 * @param th Polar angle theta.
 * @return true if static observer is timelike (outside ergosphere).
 *
 * @note Use this to validate camera placement before calling
 *       init_ray_from_camera_general(). If false, use ZAMO camera instead.
 */
__device__ __host__ inline bool kerr_static_timelike(double M, double a, double r, double th) {
    double g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov_d(M, a, r, th, g_tt, g_tphi, g_phiphi, g_rr, g_thth);
    return g_tt < 0.0;
}

/**
 * @brief Initialize ray from ZAMO (Zero Angular Momentum Observer) camera.
 *
 * ZAMOs have zero angular momentum (L = 0) and are dragged along by the
 * black hole's rotation. They exist everywhere outside the horizon, including
 * inside the ergosphere where static observers cannot exist.
 *
 * The ZAMO 4-velocity is:
 *   u^t = 1/alpha
 *   u^phi = omega/alpha
 * where alpha is the lapse and omega is the frame-dragging angular velocity.
 *
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param r_cam Camera radial position.
 * @param theta_cam Camera polar angle.
 * @param phi_cam Camera azimuthal angle.
 * @param n_r Local radial direction component.
 * @param n_theta Local polar direction component.
 * @param n_phi Local azimuthal direction component.
 * @return Initialized GeodesicStateD ready for integration.
 *
 * @note Valid inside and outside the ergosphere.
 */
__device__ __host__ inline GeodesicStateD init_ray_from_camera_zamo_d(double M, double a,
                                                                     double r_cam,
                                                                     double theta_cam,
                                                                     double phi_cam,
                                                                     double n_r, double n_theta, double n_phi) {
    GeodesicStateD s{};
    s.t = 0.0; s.r = r_cam; s.theta = theta_cam; s.phi = phi_cam;

    // Compute metric scalars in double precision
    double ct = cos(theta_cam), st = sin(theta_cam);
    double r2 = r_cam * r_cam, a2 = a * a;
    double Sigma = fmax(BH_EPS_D, r2 + a2 * ct * ct);
    double Delta = r2 - 2.0 * M * r_cam + a2;
    double A = fmax(BH_EPS_D, (r2 + a2) * (r2 + a2) - a2 * Delta * st * st);

    // Lapse and frame-dragging angular velocity (with A protection)
    double alpha = sqrt(fmax(BH_EPS_D, (Delta * Sigma) / A));
    double omega = (2.0 * M * a * r_cam) / A;

    // ZAMO tetrad basis vectors (contravariant components)
    // e_(t) = (1/alpha) (partial_t + omega partial_phi)
    double e_t_t = 1.0 / fmax(BH_EPS_D, alpha);
    double e_t_phi = omega / fmax(BH_EPS_D, alpha);

    // e_(r) = sqrt(Delta/Sigma) partial_r
    double e_r_r = sqrt(fmax(BH_EPS_D, Delta / Sigma));

    // e_(theta) = (1/sqrt(Sigma)) partial_theta
    double e_th_th = 1.0 / sqrt(Sigma);

    // e_(phi) = sqrt(Sigma/A) / sin(theta) partial_phi
    double e_phi_phi = sqrt(fmax(BH_EPS_D, Sigma / A)) / fmax(BH_EPS_D, fabs(st));

    // Construct contravariant 4-momentum from local direction.
    // For a photon with local direction (n_r, n_theta, n_phi) in the ZAMO frame:
    //   p^mu = E_loc * [e^mu_(t) + n_r * e^mu_(r) + n_theta * e^mu_(theta) + n_phi * e^mu_(phi)]
    //
    // The timelike leg e_(t) contributes to BOTH p^t and p^phi (frame dragging):
    //   p^t   = e^t_(t) = 1/alpha
    //   p^phi = e^phi_(t) + n_phi * e^phi_(phi) = omega/alpha + n_phi * sqrt(Sigma/A)/sin(theta)
    //
    // Note: Even for n_phi = 0, p^phi != 0 because the ZAMO co-rotates with spacetime.
    double pt_contra = e_t_t;
    double pr_contra = n_r * e_r_r;
    double pth_contra = n_theta * e_th_th;
    double pph_contra = e_t_phi + n_phi * e_phi_phi;

    // Lower indices using true FP64 covariant metric
    double g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov_d(M, a, r_cam, theta_cam, g_tt, g_tphi, g_phiphi, g_rr, g_thth);

    s.pt = g_tt * pt_contra + g_tphi * pph_contra;
    s.pr = g_rr * pr_contra;
    s.ptheta = g_thth * pth_contra;
    s.pphi = g_tphi * pt_contra + g_phiphi * pph_contra;

    return s;
}

/* ============================================================================
 * ADAPTIVE RK4 INTEGRATOR (FP64)
 * ============================================================================
 */

/**
 * @brief Perform a single RK4 step (helper for adaptive integration).
 *
 * Takes one RK4 step of size h and returns the result in s_out.
 * Also checks termination conditions.
 *
 * @param s Input state (not modified).
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param h Step size.
 * @param[out] s_out Output state after step.
 * @return true if integration can continue, false if terminated.
 */
__device__ inline bool rk4_step_d(const GeodesicStateD& s, double M, double a, double h, GeodesicStateD& s_out) {
    auto deriv = [&](const GeodesicStateD& st, GeodesicStateD& d) {
        computeGeodesicDerivativesKerrD(st, M, a, d);
    };

    GeodesicStateD k1, k2, k3, k4, tmp, st = s;

    deriv(st, k1);

    tmp.t = st.t + 0.5*h*k1.t; tmp.r = st.r + 0.5*h*k1.r;
    tmp.theta = st.theta + 0.5*h*k1.theta; tmp.phi = st.phi + 0.5*h*k1.phi;
    tmp.pt = st.pt + 0.5*h*k1.pt; tmp.pr = st.pr + 0.5*h*k1.pr;
    tmp.ptheta = st.ptheta + 0.5*h*k1.ptheta; tmp.pphi = st.pphi + 0.5*h*k1.pphi;

    deriv(tmp, k2);

    tmp.t = st.t + 0.5*h*k2.t; tmp.r = st.r + 0.5*h*k2.r;
    tmp.theta = st.theta + 0.5*h*k2.theta; tmp.phi = st.phi + 0.5*h*k2.phi;
    tmp.pt = st.pt + 0.5*h*k2.pt; tmp.pr = st.pr + 0.5*h*k2.pr;
    tmp.ptheta = st.ptheta + 0.5*h*k2.ptheta; tmp.pphi = st.pphi + 0.5*h*k2.pphi;

    deriv(tmp, k3);

    tmp.t = st.t + h*k3.t; tmp.r = st.r + h*k3.r;
    tmp.theta = st.theta + h*k3.theta; tmp.phi = st.phi + h*k3.phi;
    tmp.pt = st.pt + h*k3.pt; tmp.pr = st.pr + h*k3.pr;
    tmp.ptheta = st.ptheta + h*k3.ptheta; tmp.pphi = st.pphi + h*k3.pphi;

    deriv(tmp, k4);

    double h6 = h / 6.0;
    s_out = st;
    s_out.t += h6 * bh_rk4_sum(k1.t, k2.t, k3.t, k4.t);
    s_out.r += h6 * bh_rk4_sum(k1.r, k2.r, k3.r, k4.r);
    s_out.theta += h6 * bh_rk4_sum(k1.theta, k2.theta, k3.theta, k4.theta);
    s_out.phi += h6 * bh_rk4_sum(k1.phi, k2.phi, k3.phi, k4.phi);
    s_out.pt += h6 * bh_rk4_sum(k1.pt, k2.pt, k3.pt, k4.pt);
    s_out.pr += h6 * bh_rk4_sum(k1.pr, k2.pr, k3.pr, k4.pr);
    s_out.ptheta += h6 * bh_rk4_sum(k1.ptheta, k2.ptheta, k3.ptheta, k4.ptheta);
    s_out.pphi += h6 * bh_rk4_sum(k1.pphi, k2.pphi, k3.pphi, k4.pphi);

    double disc = M*M - a*a;
    if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    return (s_out.r > r_plus + 1e-12 && s_out.theta > 1e-9 && s_out.theta < BH_PI_D - 1e-9);
}

/**
 * @brief Perform one adaptive RK4 step using step-doubling error estimation.
 *
 * Estimates the local truncation error by comparing:
 * - One step of size h
 * - Two steps of size h/2
 *
 * If the error is within tolerance, accepts the step and suggests a factor
 * to increase h for the next step. If the error is too large, shrinks h
 * and returns 0 to indicate rejection.
 *
 * Error tolerance: |y_full - y_half| <= abs_tol + rel_tol * max(|y_full|, |y_half|)
 *
 * @param[in,out] s Geodesic state. Updated to the new state if step accepted.
 * @param M Black hole mass.
 * @param a Spin parameter.
 * @param[in,out] h Step size. May be reduced on rejection.
 * @param rel_tol Relative tolerance (default: BH_RK_REL_TOL = 1e-10).
 * @param abs_tol Absolute tolerance (default: BH_RK_ABS_TOL = 1e-12).
 * @return Factor to multiply h for next step (>0 if accepted, 0 if rejected).
 *
 * @note On acceptance, s is updated to the more accurate half-step result.
 * @note On rejection, s is unchanged but h is reduced.
 */
__device__ inline double rk4_adaptive_step_d(GeodesicStateD& s, double M, double a, double& h,
                                             double rel_tol = BH_RK_REL_TOL, double abs_tol = BH_RK_ABS_TOL) {
    // Compute full step and two half-steps
    GeodesicStateD s_full, s_half1, s_half2;
    bool ok1 = rk4_step_d(s, M, a, h, s_full);
    double h2 = 0.5 * h;
    bool ok2 = rk4_step_d(s, M, a, h2, s_half1);
    bool ok3 = rk4_step_d(s_half1, M, a, h2, s_half2);

    // If any step failed (hit horizon/pole), shrink and signal rejection
    if (!(ok1 && ok2 && ok3)) {
        h *= 0.5;
        return 0.0;  // Reject step - state s is unchanged
    }

    // Compute normalized error for each component
    auto err = [](double a, double b, double scale) {
        return fabs(a - b) / fmax(scale, 1e-30);
    };

    double scale_r = abs_tol + rel_tol * fmax(fabs(s_full.r), fabs(s_half2.r));
    double scale_th = abs_tol + rel_tol * fmax(fabs(s_full.theta), fabs(s_half2.theta));
    double scale_ph = abs_tol + rel_tol * fmax(fabs(s_full.phi), fabs(s_half2.phi));
    double scale_pr = abs_tol + rel_tol * fmax(fabs(s_full.pr), fabs(s_half2.pr));
    double scale_pth = abs_tol + rel_tol * fmax(fabs(s_full.ptheta), fabs(s_half2.ptheta));
    double scale_pph = abs_tol + rel_tol * fmax(fabs(s_full.pphi), fabs(s_half2.pphi));

    // Maximum normalized error across all components
    double e = fmax(
        fmax(err(s_full.r, s_half2.r, scale_r), err(s_full.theta, s_half2.theta, scale_th)),
        fmax(
            fmax(err(s_full.phi, s_half2.phi, scale_ph), err(s_full.pr, s_half2.pr, scale_pr)),
            fmax(err(s_full.ptheta, s_half2.ptheta, scale_pth), err(s_full.pphi, s_half2.pphi, scale_pph))
        )
    );

    // Accept step if error within tolerance
    if (e <= 1.0) {
        s = s_half2;  // Use more accurate result
        // Suggest step increase for next iteration (capped at 2x)
        return fmin(2.0, BH_RK_SAFETY * pow(fmax(1e-12, e), -0.2));
    }

    // Reject step - shrink h and signal rejection
    double fac = fmax(0.1, BH_RK_SAFETY * pow(fmax(1e-12, e), -0.25));
    h *= fac;
    return 0.0;
}

#pragma endregion

#endif // BH_KERR_H
