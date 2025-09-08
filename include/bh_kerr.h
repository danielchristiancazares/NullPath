// Header-only Kerr metric blocks, derivatives, Hamiltonian, RK4, and camera init
#ifndef BH_KERR_H
#define BH_KERR_H

#include "bh_common.h"
#include "bh_types.h"

// Kerr metric helpers (Boyer-Lindquist), covariant and contravariant blocks
__device__ __host__ inline void kerr_blocks_cov(float M, float a, float r, float th,
                                               float& g_tt, float& g_tphi, float& g_phiphi,
                                               float& g_rr, float& g_thth) {
    float ct = cosf(th);
    float st = sinf(th);
    float st2 = st*st;
    float r2 = r*r;
    float a2 = a*a;
    float Sigma = r2 + a2 * ct*ct;
    float Delta = r2 - 2.0f*M*r + a2;
    float A = (r2 + a2)*(r2 + a2) - a2 * Delta * st2;
    g_tt    = -(1.0f - 2.0f*M*r / Sigma);
    g_tphi  = -2.0f * M * a * r * st2 / Sigma;
    g_phiphi=  A * st2 / Sigma;
    g_rr    = Sigma / fmaxf(BH_EPS, Delta);
    g_thth  = Sigma;
}

__device__ __host__ inline void kerr_blocks_contra(float M, float a, float r, float th,
                                                  float& gtt, float& gtphi, float& gphiphi,
                                                  float& grr, float& gthth) {
    float ct = cosf(th);
    float st = sinf(th);
    float st2 = st*st;
    float r2 = r*r;
    float a2 = a*a;
    float Sigma = r2 + a2 * ct*ct;
    float Delta = r2 - 2.0f*M*r + a2;
    float A = (r2 + a2)*(r2 + a2) - a2 * Delta * st2;
    grr   = fmaxf(BH_EPS, Delta) / fmaxf(BH_EPS, Sigma);
    gthth = 1.0f / fmaxf(BH_EPS, Sigma);
    gtt   = -A / fmaxf(BH_EPS, Sigma * Delta);
    gtphi = -2.0f * M * a * r / fmaxf(BH_EPS, Sigma * Delta);
    gphiphi = (Delta - a2 * st2) / fmaxf(BH_EPS, Sigma * Delta * st2 + BH_EPS);
}

// Analytic derivatives of inverse metric components wrt r and theta
__device__ __host__ inline void kerr_blocks_contra_derivs(float M, float a, float r, float th,
                                                        float& gtt_r, float& gtphi_r, float& gphiphi_r,
                                                        float& grr_r, float& gthth_r,
                                                        float& gtt_th, float& gtphi_th, float& gphiphi_th,
                                                        float& grr_th, float& gthth_th) {
    float s = sinf(th), c = cosf(th);
    float s2 = fmaxf(BH_EPS_SIN, s*s);
    float r2 = r*r;
    float a2 = a*a;
    float B = r2 + a2;
    float Sigma = r2 + a2 * c*c;
    float Delta = r2 - 2.0f*M*r + a2;
    float A = (B*B) - a2 * Delta * s2;
    float D = Sigma * Delta;
    float D2 = fmaxf(1e-24f, D*D);
    float s2_theta = 2.0f * s * c;
    // Partials of Sigma, Delta, A, D
    float Sigma_r = 2.0f * r;
    float Sigma_th = -2.0f * a2 * s * c;
    float Delta_r = 2.0f * (r - M);
    float A_r = 4.0f * r * B - a2 * Delta_r * s2; // (d/dr)(B^2) = 2B*2r
    float A_th = - a2 * Delta * s2_theta;
    float D_r = Sigma_r * Delta + Sigma * Delta_r;
    float D_th = Sigma_th * Delta; // Delta_th = 0

    // g^{tt} = -A/D
    gtt_r = - (A_r * D - A * D_r) / D2;
    gtt_th = - (A_th * D - A * D_th) / D2;

    // g^{tφ} = - 2 a M r / D
    float N_tphi = 2.0f * a * M * r;
    float N_tphi_r = 2.0f * a * M;
    float N_tphi_th = 0.0f;
    gtphi_r = - (N_tphi_r * D - N_tphi * D_r) / D2;
    gtphi_th = - (N_tphi_th * D - N_tphi * D_th) / D2; // simplifies to +N*D_th/D^2

    // g^{rr} = Delta / Sigma
    grr_r = (Delta_r * Sigma - Delta * Sigma_r) / (Sigma * Sigma);
    grr_th = - (Delta * Sigma_th) / (Sigma * Sigma);

    // g^{θθ} = 1 / Sigma
    gthth_r = - Sigma_r / (Sigma * Sigma);
    gthth_th = - Sigma_th / (Sigma * Sigma);

    // g^{φφ} = (Delta - a^2 s^2) / (Delta Sigma s^2)
    float N = Delta - a2 * s2;
    float N_r = Delta_r;
    float N_th = - a2 * s2_theta;
    float Den = Delta * Sigma * s2;
    float Den_r = Delta_r * Sigma * s2 + Delta * Sigma_r * s2; // s2 independent of r
    float Den_th = Delta * Sigma_th * s2 + Delta * Sigma * s2_theta;
    float Den2 = fmaxf(1e-24f, Den*Den);
    gphiphi_r = (N_r * Den - N * Den_r) / Den2;
    gphiphi_th = (N_th * Den - N * Den_th) / Den2;
}

// Hamiltonian for Kerr null geodesics (contravariant form)
__host__ __device__ inline float hamiltonian_kerr(const GeodesicState& s, float M, float a) {
    float gtt, gtphi, gphiphi, grr, gthth;
    kerr_blocks_contra(M, a, s.r, s.theta, gtt, gtphi, gphiphi, grr, gthth);
    float cross = 2.0f * gtphi * s.pt * s.pphi;
    float H = 0.5f * (gtt*s.pt*s.pt + cross + gphiphi*s.pphi*s.pphi + grr*s.pr*s.pr + gthth*s.ptheta*s.ptheta);
    return H;
}

// Kerr geodesic derivatives using Hamiltonian form and analytic partials
__device__ inline void computeGeodesicDerivativesKerr(const GeodesicState& st, float M, float a,
                                                      GeodesicState& d) {
    float r = st.r;
    float th = st.theta;
    // Horizon radius r_+ = M + sqrt(M^2 - a^2)
    float disc = fmaxf(0.0f, M*M - a*a);
    float r_plus = M + sqrtf(disc);
    if (r <= r_plus + 1e-3f || r <= 1e-6f) { d = GeodesicState(); return; }

    float gtt, gtphi, gphiphi, grr, gthth;
    kerr_blocks_contra(M, a, r, th, gtt, gtphi, gphiphi, grr, gthth);
    float gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r;
    float gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th;
    kerr_blocks_contra_derivs(M, a, r, th,
                              gtt_r, gtphi_r, gphiphi_r, grr_r, gthth_r,
                              gtt_th, gtphi_th, gphiphi_th, grr_th, gthth_th);

    // dx/dλ = g^{μν} p_ν
    d.t     = gtt   * st.pt + gtphi   * st.pphi;
    d.r     = grr   * st.pr;
    d.theta = gthth * st.ptheta;
    d.phi   = gtphi * st.pt + gphiphi * st.pphi;

    // dp/dλ = -1/2 ∂ g^{αβ} p_α p_β
    d.pt = 0.0f;   // stationarity
    d.pphi = 0.0f; // axisymmetry
    float S_r = gtt_r   * st.pt*st.pt + 2.0f*gtphi_r * st.pt*st.pphi + gphiphi_r * st.pphi*st.pphi
                + grr_r * st.pr*st.pr + gthth_r * st.ptheta*st.ptheta;
    float S_th= gtt_th  * st.pt*st.pt + 2.0f*gtphi_th* st.pt*st.pphi + gphiphi_th* st.pphi*st.pphi
                + grr_th* st.pr*st.pr + gthth_th* st.ptheta*st.ptheta;
    d.pr = -0.5f * S_r;
    d.ptheta = -0.5f * S_th;
}

__device__ inline bool integrateGeodesicKerr(GeodesicState& st, float M, float a, float h) {
    GeodesicState k1, k2, k3, k4, tmp;
    computeGeodesicDerivativesKerr(st, M, a, k1);
    tmp.t = st.t + 0.5f*h*k1.t; tmp.r = st.r + 0.5f*h*k1.r; tmp.theta = st.theta + 0.5f*h*k1.theta; tmp.phi = st.phi + 0.5f*h*k1.phi;
    tmp.pt = st.pt + 0.5f*h*k1.pt; tmp.pr = st.pr + 0.5f*h*k1.pr; tmp.ptheta = st.ptheta + 0.5f*h*k1.ptheta; tmp.pphi = st.pphi + 0.5f*h*k1.pphi;
    computeGeodesicDerivativesKerr(tmp, M, a, k2);
    tmp.t = st.t + 0.5f*h*k2.t; tmp.r = st.r + 0.5f*h*k2.r; tmp.theta = st.theta + 0.5f*h*k2.theta; tmp.phi = st.phi + 0.5f*h*k2.phi;
    tmp.pt = st.pt + 0.5f*h*k2.pt; tmp.pr = st.pr + 0.5f*h*k2.pr; tmp.ptheta = st.ptheta + 0.5f*h*k2.ptheta; tmp.pphi = st.pphi + 0.5f*h*k2.pphi;
    computeGeodesicDerivativesKerr(tmp, M, a, k3);
    tmp.t = st.t + h*k3.t; tmp.r = st.r + h*k3.r; tmp.theta = st.theta + h*k3.theta; tmp.phi = st.phi + h*k3.phi;
    tmp.pt = st.pt + h*k3.pt; tmp.pr = st.pr + h*k3.pr; tmp.ptheta = st.ptheta + h*k3.ptheta; tmp.pphi = st.pphi + h*k3.pphi;
    computeGeodesicDerivativesKerr(tmp, M, a, k4);
    float h6 = h/6.0f;
    st.t += h6 * (k1.t + 2*k2.t + 2*k3.t + k4.t);
    st.r += h6 * (k1.r + 2*k2.r + 2*k3.r + k4.r);
    st.theta += h6 * (k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
    st.phi += h6 * (k1.phi + 2*k2.phi + 2*k3.phi + k4.phi);
    st.pt += h6 * (k1.pt + 2*k2.pt + 2*k3.pt + k4.pt);
    st.pr += h6 * (k1.pr + 2*k2.pr + 2*k3.pr + k4.pr);
    st.ptheta += h6 * (k1.ptheta + 2*k2.ptheta + 2*k3.ptheta + k4.ptheta);
    st.pphi += h6 * (k1.pphi + 2*k2.pphi + 2*k3.pphi + k4.pphi);
    // bounds
    float disc = fmaxf(0.0f, M*M - a*a);
    float r_plus = M + sqrtf(disc);
    return (st.r > r_plus + 1e-3f && st.theta > 1e-3f && st.theta < (float)BH_PI - 1e-3f);
}

// Camera ray init for static observer (valid outside ergosphere)
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
        float norm2 = g_tt + 2.0f * g_tphi * beta + g_phiphi * beta * beta;
        float inv = rsqrtf(fabsf(norm2) + 1e-12f);
        ephi_t = inv; ephi_phi = beta * inv;
    } else {
        ephi_t = 0.0f; ephi_phi = rsqrtf(fmaxf(BH_EPS, g_phiphi));
    }

    const float E_loc = 1.0f;
    float pt_contra = E_loc * (ut + n_phi * ephi_t);
    float pr_contra = E_loc * (n_r * er_r);
    float pth_contra = E_loc * (n_theta * eth_th);
    float pph_contra = E_loc * (n_phi * ephi_phi);

    s.pt = g_tt * pt_contra + g_tphi * pph_contra;
    s.pr = g_rr * pr_contra;
    s.ptheta = g_thth * pth_contra;
    s.pphi = g_tphi * pt_contra + g_phiphi * pph_contra;
    return s;
}

// Emitter 4-velocity (circular, equatorial) and Doppler g-factor
__device__ inline float doppler_g_general(const GeodesicState& st, float M, float a) {
    float r = st.r;
    float th = (float)BH_PI * 0.5f; // disk plane
    float g_tt, g_tphi, g_phiphi, g_rr, g_thth;
    kerr_blocks_cov(M, a, r, th, g_tt, g_tphi, g_phiphi, g_rr, g_thth);
    // Prograde Keplerian angular velocity (BL coords)
    float r_over_M = r / fmaxf(BH_EPS, M);
    float Omega = 1.0f / (a / fmaxf(BH_EPS, M) + powf(fmaxf(1.0f, r_over_M), 1.5f)) / fmaxf(BH_EPS, M);
    float denom = -(g_tt + 2.0f * g_tphi * Omega + g_phiphi * Omega * Omega);
    float u_t = rsqrtf(fmaxf(BH_EPS, denom));
    float u_phi = Omega * u_t;
    float E_inf = -st.pt;
    float k_u = -(st.pt * u_t + st.pphi * u_phi);
    float g = E_inf / fmaxf(1e-9f, k_u);
    return fmaxf(0.0f, g);
}

#pragma region FP64_geodesics
// ===================== FP64 (double) geodesic math =====================
struct GeodesicStateD { double t,r,theta,phi, pt,pr,ptheta,pphi; };

__device__ __host__ inline void kerr_blocks_contra_d(double M, double a, double r, double th,
                                                    double& gtt, double& gtphi, double& gphiphi,
                                                    double& grr, double& gthth) {
    double ct = cos(th);
    double st = sin(th);
    double st2 = st*st > 1e-30 ? st*st : 1e-30;
    double r2 = r*r;
    double a2 = a*a;
    double Sigma = r2 + a2 * ct*ct;
    double Delta = r2 - 2.0*M*r + a2;
    double A = (r2 + a2)*(r2 + a2) - a2 * Delta * st2;
    grr   = Delta / Sigma;
    gthth = 1.0 / Sigma;
    gtt   = -A / (Sigma * Delta);
    gtphi = -2.0 * M * a * r / (Sigma * Delta);
    gphiphi = (Delta - a2 * st2) / (Sigma * Delta * st2);
}

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
    double A_th = - a2 * Delta * s2_th;
    double D_r = Sigma_r * Delta + Sigma * Delta_r;
    double D_th = Sigma_th * Delta;
    gtt_r = - (A_r * D - A * D_r) / D2;
    gtt_th = - (A_th * D - A * D_th) / D2;
    double Nt = 2.0 * M * a * r;
    double Nt_r = 2.0 * M * a;
    gtphi_r = - (Nt_r * D - Nt * D_r) / D2;
    gtphi_th = - (0.0 * D - Nt * D_th) / D2;
    grr_r = (Delta_r * Sigma - Delta * Sigma_r) / (Sigma * Sigma);
    grr_th = - (Delta * Sigma_th) / (Sigma * Sigma);
    gthth_r = - Sigma_r / (Sigma * Sigma);
    gthth_th = - Sigma_th / (Sigma * Sigma);
    double N = Delta - a2 * s2;
    double N_r = Delta_r;
    double N_th = - a2 * s2_th;
    double Den = Delta * Sigma * s2 + 1e-48;
    double Den_r = Delta_r * Sigma * s2 + Delta * Sigma_r * s2;
    double Den_th = Delta * Sigma_th * s2 + Delta * Sigma * s2_th;
    double Den2 = Den * Den + 1e-48;
    gphiphi_r = (N_r * Den - N * Den_r) / Den2;
    gphiphi_th = (N_th * Den - N * Den_th) / Den2;
}

__device__ inline void computeGeodesicDerivativesKerrD(const GeodesicStateD& st, double M, double a,
                                                       GeodesicStateD& d) {
    double r = st.r;
    double th = st.theta;
    double disc = M*M - a*a; if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    if (r <= r_plus + 1e-12 || r <= 1e-12) { d = GeodesicStateD(); return; }
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
    d.pt = 0.0; d.pphi = 0.0;
    double S_r = gtt_r   * st.pt*st.pt + 2.0*gtphi_r * st.pt*st.pphi + gphiphi_r * st.pphi*st.pphi
                 + grr_r * st.pr*st.pr + gthth_r * st.ptheta*st.ptheta;
    double S_th= gtt_th  * st.pt*st.pt + 2.0*gtphi_th* st.pt*st.pphi + gphiphi_th* st.pphi*st.pphi
                 + grr_th* st.pr*st.pr + gthth_th* st.ptheta*st.ptheta;
    d.pr = -0.5 * S_r;
    d.ptheta = -0.5 * S_th;
}

__device__ inline bool integrateGeodesicKerrD(GeodesicStateD& st, double M, double a, double h) {
    GeodesicStateD k1, k2, k3, k4, tmp;
    computeGeodesicDerivativesKerrD(st, M, a, k1);
    tmp.t = st.t + 0.5*h*k1.t; tmp.r = st.r + 0.5*h*k1.r; tmp.theta = st.theta + 0.5*h*k1.theta; tmp.phi = st.phi + 0.5*h*k1.phi;
    tmp.pt = st.pt + 0.5*h*k1.pt; tmp.pr = st.pr + 0.5*h*k1.pr; tmp.ptheta = st.ptheta + 0.5*h*k1.ptheta; tmp.pphi = st.pphi + 0.5*h*k1.pphi;
    computeGeodesicDerivativesKerrD(tmp, M, a, k2);
    tmp.t = st.t + 0.5*h*k2.t; tmp.r = st.r + 0.5*h*k2.r; tmp.theta = st.theta + 0.5*h*k2.theta; tmp.phi = st.phi + 0.5*h*k2.phi;
    tmp.pt = st.pt + 0.5*h*k2.pt; tmp.pr = st.pr + 0.5*h*k2.pr; tmp.ptheta = st.ptheta + 0.5*h*k2.ptheta; tmp.pphi = st.pphi + 0.5*h*k2.pphi;
    computeGeodesicDerivativesKerrD(tmp, M, a, k3);
    tmp.t = st.t + h*k3.t; tmp.r = st.r + h*k3.r; tmp.theta = st.theta + h*k3.theta; tmp.phi = st.phi + h*k3.phi;
    tmp.pt = st.pt + h*k3.pt; tmp.pr = st.pr + h*k3.pr; tmp.ptheta = st.ptheta + h*k3.ptheta; tmp.pphi = st.pphi + h*k3.pphi;
    computeGeodesicDerivativesKerrD(tmp, M, a, k4);
    double h6 = h/6.0;
    st.t += h6 * (k1.t + 2*k2.t + 2*k3.t + k4.t);
    st.r += h6 * (k1.r + 2*k2.r + 2*k3.r + k4.r);
    st.theta += h6 * (k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
    st.phi += h6 * (k1.phi + 2*k2.phi + 2*k3.phi + k4.phi);
    st.pt += h6 * (k1.pt + 2*k2.pt + 2*k3.pt + k4.pt);
    st.pr += h6 * (k1.pr + 2*k2.pr + 2*k3.pr + k4.pr);
    st.ptheta += h6 * (k1.ptheta + 2*k2.ptheta + 2*k3.ptheta + k4.ptheta);
    st.pphi += h6 * (k1.pphi + 2*k2.pphi + 2*k3.pphi + k4.pphi);
    double disc = M*M - a*a; if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    return (st.r > r_plus + 1e-9 && st.theta > 1e-9 && st.theta < BH_PI_D - 1e-9);
}

__device__ __host__ inline GeodesicStateD init_ray_from_camera_general_d(double M, double a,
                                                                        double r_cam,
                                                                        double theta_cam,
                                                                        double phi_cam,
                                                                        double n_r, double n_theta, double n_phi) {
    GeodesicStateD s{};
    s.t = 0.0; s.r = r_cam; s.theta = theta_cam; s.phi = phi_cam;
    float g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f;
    kerr_blocks_cov((float)M, (float)a, (float)s.r, (float)s.theta, g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f);
    double g_tt = (double)g_tt_f, g_tphi = (double)g_tphi_f, g_phiphi = (double)g_phiphi_f, g_rr = (double)g_rr_f, g_thth = (double)g_thth_f;
    double ut = rsqrt(-(g_tt) > 1e-24 ? -(g_tt) : 1e-24);
    double er_r = rsqrt(g_rr > 1e-24 ? g_rr : 1e-24);
    double eth_th = rsqrt(g_thth > 1e-24 ? g_thth : 1e-24);
    double beta = (fabs(g_tphi) > 1e-24) ? (-g_tt / g_tphi) : 0.0;
    double norm2 = g_tt + 2.0 * g_tphi * beta + g_phiphi * beta * beta;
    double inv = rsqrt(fabs(norm2) > 1e-24 ? fabs(norm2) : 1e-24);
    double ephi_t = inv, ephi_phi = beta * inv;
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

__device__ inline double doppler_g_general_d(const GeodesicStateD& st, double M, double a) {
    double r = st.r;
    double th = BH_PI_D * 0.5; // disk plane
    float g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f;
    kerr_blocks_cov((float)M, (float)a, (float)r, (float)th, g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f);
    double g_tt = (double)g_tt_f, g_tphi = (double)g_tphi_f, g_phiphi = (double)g_phiphi_f;
    double r_over_M = r / (M > 1e-12 ? M : 1e-12);
    double Omega = 1.0 / ((a / (M > 1e-12 ? M : 1e-12)) + pow(r_over_M > 1.0 ? r_over_M : 1.0, 1.5)) / (M > 1e-12 ? M : 1.0);
    double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
    double u_t = rsqrt(denom > 1e-24 ? denom : 1e-24);
    double u_phi = Omega * u_t;
    double E_inf = -st.pt;
    double k_u = -(st.pt * u_t + st.pphi * u_phi);
    double g = E_inf / (k_u > 1e-18 ? k_u : 1e-18);
    return g > 0.0 ? g : 0.0;
}

// Static timelike check (BL): returns true if static observer (∂_t) is timelike
__device__ __host__ inline bool kerr_static_timelike(double M, double a, double r, double th) {
    float g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f;
    kerr_blocks_cov((float)M, (float)a, (float)r, (float)th, g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f);
    return g_tt_f < 0.0f;
}

// ZAMO-based camera ray init (BL coordinates)
__device__ __host__ inline GeodesicStateD init_ray_from_camera_zamo_d(double M, double a,
                                                                     double r_cam,
                                                                     double theta_cam,
                                                                     double phi_cam,
                                                                     double n_r, double n_theta, double n_phi) {
    GeodesicStateD s{}; s.t=0.0; s.r=r_cam; s.theta=theta_cam; s.phi=phi_cam;
    // ZAMO tetrad components in BL
    float Sigma_f, Delta_f, A_f, s_f, c_f;
    // Recompute scalars in double
    double ct = cos(theta_cam), st = sin(theta_cam);
    double r2 = r_cam * r_cam, a2 = a * a;
    double Sigma = r2 + a2 * ct*ct;
    double Delta = r2 - 2.0 * M * r_cam + a2;
    double A = (r2 + a2)*(r2 + a2) - a2 * Delta * st*st;
    double alpha = sqrt(fmax(1e-24, (Delta * Sigma) / A));
    double omega = (2.0 * M * a * r_cam) / A;
    // e_(t) = 1/alpha (∂_t + omega ∂_phi)
    double e_t_t = 1.0 / alpha;
    double e_t_phi = omega / alpha;
    // e_(r) = sqrt(Delta/Sigma) ∂_r ; e_(theta)= 1/sqrt(Sigma) ∂_theta ; e_(phi) = sqrt(Sigma/A) 1/sinθ ∂_phi
    double e_r_r = sqrt(fmax(1e-24, Delta / Sigma));
    double e_th_th = 1.0 / sqrt(fmax(1e-24, Sigma));
    double e_phi_phi = sqrt(fmax(1e-24, Sigma / A)) / fmax(1e-12, st);
    // Local photon direction components (unit): along -e_r + n_phi e_phi + n_theta e_theta
    double pt_contra = (1.0) * (e_t_t + n_phi * e_t_phi);
    double pr_contra = (1.0) * (n_r * e_r_r);
    double pth_contra = (1.0) * (n_theta * e_th_th);
    double pph_contra = (1.0) * (n_phi * e_phi_phi);
    // Lower with g_{μν}
    float g_tt_f, g_tphi_f, g_phiphi_f, g_rr_f, g_thth_f;
    kerr_blocks_cov((float)M,(float)a,(float)r_cam,(float)theta_cam,g_tt_f,g_tphi_f,g_phiphi_f,g_rr_f,g_thth_f);
    double g_tt = (double)g_tt_f, g_tphi = (double)g_tphi_f, g_phiphi = (double)g_phiphi_f, g_rr = (double)g_rr_f, g_thth = (double)g_thth_f;
    s.pt = g_tt * pt_contra + g_tphi * pph_contra;
    s.pr = g_rr * pr_contra;
    s.ptheta = g_thth * pth_contra;
    s.pphi = g_tphi * pt_contra + g_phiphi * pph_contra;
    return s;
}

// Adaptive RK4 by step-doubling (embedded error via one step vs two half-steps), double precision
__device__ inline bool rk4_step_d(const GeodesicStateD& s, double M, double a, double h, GeodesicStateD& s_out) {
    auto deriv = [&](const GeodesicStateD& st, GeodesicStateD& d){ computeGeodesicDerivativesKerrD(st,M,a,d); };
    GeodesicStateD k1,k2,k3,k4,tmp, st=s;
    deriv(st,k1);
    tmp.t=st.t+0.5*h*k1.t; tmp.r=st.r+0.5*h*k1.r; tmp.theta=st.theta+0.5*h*k1.theta; tmp.phi=st.phi+0.5*h*k1.phi;
    tmp.pt=st.pt+0.5*h*k1.pt; tmp.pr=st.pr+0.5*h*k1.pr; tmp.ptheta=st.ptheta+0.5*h*k1.ptheta; tmp.pphi=st.pphi+0.5*h*k1.pphi;
    deriv(tmp,k2);
    tmp.t=st.t+0.5*h*k2.t; tmp.r=st.r+0.5*h*k2.r; tmp.theta=st.theta+0.5*h*k2.theta; tmp.phi=st.phi+0.5*h*k2.phi;
    tmp.pt=st.pt+0.5*h*k2.pt; tmp.pr=st.pr+0.5*h*k2.pr; tmp.ptheta=st.ptheta+0.5*h*k2.ptheta; tmp.pphi=st.pphi+0.5*h*k2.pphi;
    deriv(tmp,k3);
    tmp.t=st.t+h*k3.t; tmp.r=st.r+h*k3.r; tmp.theta=st.theta+h*k3.theta; tmp.phi=st.phi+h*k3.phi;
    tmp.pt=st.pt+h*k3.pt; tmp.pr=st.pr+h*k3.pr; tmp.ptheta=st.ptheta+h*k3.ptheta; tmp.pphi=st.pphi+h*k3.pphi;
    deriv(tmp,k4);
    double h6 = h/6.0;
    s_out = st;
    s_out.t     += h6*(k1.t + 2*k2.t + 2*k3.t + k4.t);
    s_out.r     += h6*(k1.r + 2*k2.r + 2*k3.r + k4.r);
    s_out.theta += h6*(k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
    s_out.phi   += h6*(k1.phi + 2*k2.phi + 2*k3.phi + k4.phi);
    s_out.pt    += h6*(k1.pt + 2*k2.pt + 2*k3.pt + k4.pt);
    s_out.pr    += h6*(k1.pr + 2*k2.pr + 2*k3.pr + k4.pr);
    s_out.ptheta+= h6*(k1.ptheta + 2*k2.ptheta + 2*k3.ptheta + k4.ptheta);
    s_out.pphi  += h6*(k1.pphi + 2*k2.pphi + 2*k3.pphi + k4.pphi);
    // bounds check
    double disc = M*M - a*a; if (disc < 0.0) disc = 0.0;
    double r_plus = M + sqrt(disc);
    return (s_out.r > r_plus + 1e-12 && s_out.theta > 1e-9 && s_out.theta < BH_PI_D - 1e-9);
}

__device__ inline double rk4_adaptive_step_d(GeodesicStateD& s, double M, double a, double& h,
                                             double rel_tol = BH_RK_REL_TOL, double abs_tol = BH_RK_ABS_TOL) {
    // Try one full step h and two half-steps h/2; compare to estimate error
    GeodesicStateD s_full, s_half1, s_half2;
    bool ok1 = rk4_step_d(s, M, a, h, s_full);
    double h2 = 0.5*h;
    bool ok2 = rk4_step_d(s, M, a, h2, s_half1);
    bool ok3 = rk4_step_d(s_half1, M, a, h2, s_half2);
    if (!(ok1 && ok2 && ok3)) { h *= 0.5; return 1.0; }
    // Error norm over selected components
    auto err = [&](double a, double b, double scale){ return fabs(a-b) / fmax(scale, 1e-30); };
    double scale_r = abs_tol + rel_tol * fmax(fabs(s_full.r), fabs(s_half2.r));
    double scale_th= abs_tol + rel_tol * fmax(fabs(s_full.theta), fabs(s_half2.theta));
    double scale_ph= abs_tol + rel_tol * fmax(fabs(s_full.phi), fabs(s_half2.phi));
    double scale_pr= abs_tol + rel_tol * fmax(fabs(s_full.pr), fabs(s_half2.pr));
    double scale_pth=abs_tol + rel_tol * fmax(fabs(s_full.ptheta), fabs(s_half2.ptheta));
    double scale_pph=abs_tol + rel_tol * fmax(fabs(s_full.pphi), fabs(s_half2.pphi));
    double e = fmax(fmax(err(s_full.r, s_half2.r, scale_r), err(s_full.theta, s_half2.theta, scale_th)),
                    fmax(fmax(err(s_full.phi, s_half2.phi, scale_ph), err(s_full.pr, s_half2.pr, scale_pr)),
                         fmax(err(s_full.ptheta, s_half2.ptheta, scale_pth), err(s_full.pphi, s_half2.pphi, scale_pph))));
    // Accept if e <= 1; return suggested factor
    if (e <= 1.0) { s = s_half2; return fmin(2.0, BH_RK_SAFETY * pow(fmax(1e-12, e), -0.2)); }
    // Reject; shrink
    double fac = fmax(0.1, BH_RK_SAFETY * pow(fmax(1e-12, e), -0.25));
    h *= fac;
    return 0.0;
}

#pragma endregion

#endif // BH_KERR_H
