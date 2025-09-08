// Numerical self-test for analytic Kerr contravariant metric derivatives
// Compares kerr_blocks_contra_derivs (analytic) vs. 5-point finite differences
// across random (r, theta) samples for spins a/M in {0, 0.5, 0.9}.
// Host-only test: uses __host__ paths from bh_kerr.h (no GPU required).

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

#include "../include/bh_common.h"
#include "../include/bh_kerr.h"

struct RelErrStats {
    double max_rel_gtt_r = 0, max_rel_gtphi_r = 0, max_rel_gphiphi_r = 0, max_rel_grr_r = 0, max_rel_gthth_r = 0;
    double max_rel_gtt_th = 0, max_rel_gtphi_th = 0, max_rel_gphiphi_th = 0, max_rel_grr_th = 0, max_rel_gthth_th = 0;
    int fails = 0;
};

static inline double rel_err(double a, double b) {
    double den = std::max(1e-20, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) / den;
}

static void fd_derivs(double M, double a, double r, double th, double hr, double hth,
                      double &gtt_r, double &gtphi_r, double &gphiphi_r, double &grr_r, double &gthth_r,
                      double &gtt_th, double &gtphi_th, double &gphiphi_th, double &grr_th, double &gthth_th) {
    // r-derivatives via 5-point stencil
    auto eval = [&](double rr, double tt, double &gtt, double &gtphi, double &gphiphi, double &grr, double &gthth){
        float gtt_f, gtphi_f, gphiphi_f, grr_f, gthth_f;
        kerr_blocks_contra((float)M, (float)a, (float)rr, (float)tt, gtt_f, gtphi_f, gphiphi_f, grr_f, gthth_f);
        gtt = gtt_f; gtphi = gtphi_f; gphiphi = gphiphi_f; grr = grr_f; gthth = gthth_f;
    };
    double gtt_p2, gtphi_p2, gphiphi_p2, grr_p2, gthth_p2;
    double gtt_p1, gtphi_p1, gphiphi_p1, grr_p1, gthth_p1;
    double gtt_m1, gtphi_m1, gphiphi_m1, grr_m1, gthth_m1;
    double gtt_m2, gtphi_m2, gphiphi_m2, grr_m2, gthth_m2;

    eval(r+2*hr, th, gtt_p2, gtphi_p2, gphiphi_p2, grr_p2, gthth_p2);
    eval(r+1*hr, th, gtt_p1, gtphi_p1, gphiphi_p1, grr_p1, gthth_p1);
    eval(r-1*hr, th, gtt_m1, gtphi_m1, gphiphi_m1, grr_m1, gthth_m1);
    eval(r-2*hr, th, gtt_m2, gtphi_m2, gphiphi_m2, grr_m2, gthth_m2);

    auto fd = [&](double f_p2, double f_p1, double f_m1, double f_m2){
        return (-f_p2 + 8.0*f_p1 - 8.0*f_m1 + f_m2) / (12.0*hr);
    };
    gtt_r     = fd(gtt_p2,    gtt_p1,    gtt_m1,    gtt_m2);
    gtphi_r   = fd(gtphi_p2,  gtphi_p1,  gtphi_m1,  gtphi_m2);
    gphiphi_r = fd(gphiphi_p2,gphiphi_p1,gphiphi_m1,gphiphi_m2);
    grr_r     = fd(grr_p2,    grr_p1,    grr_m1,    grr_m2);
    gthth_r   = fd(gthth_p2,  gthth_p1,  gthth_m1,  gthth_m2);

    // theta-derivatives via 5-point stencil
    eval(r, th+2*hth, gtt_p2, gtphi_p2, gphiphi_p2, grr_p2, gthth_p2);
    eval(r, th+1*hth, gtt_p1, gtphi_p1, gphiphi_p1, grr_p1, gthth_p1);
    eval(r, th-1*hth, gtt_m1, gtphi_m1, gphiphi_m1, grr_m1, gthth_m1);
    eval(r, th-2*hth, gtt_m2, gtphi_m2, gphiphi_m2, grr_m2, gthth_m2);

    gtt_th     = fd(gtt_p2,    gtt_p1,    gtt_m1,    gtt_m2);
    gtphi_th   = fd(gtphi_p2,  gtphi_p1,  gtphi_m1,  gtphi_m2);
    gphiphi_th = fd(gphiphi_p2,gphiphi_p1,gphiphi_m1,gphiphi_m2);
    grr_th     = fd(grr_p2,    grr_p1,    grr_m1,    grr_m2);
    gthth_th   = fd(gthth_p2,  gthth_p1,  gthth_m1,  gthth_m2);
}

static RelErrStats run_suite(double M, double a, int samples, double tol, double abs_tol, unsigned seed) {
    std::mt19937 rng(seed);
    RelErrStats st;
    // Domain: r \in [r_plus + margin, 50 M], theta \in [0.05\pi, 0.95\pi]
    double disc = std::max(0.0, M*M - a*a);
    double r_plus = M + std::sqrt(disc);
    std::uniform_real_distribution<double> ur(0.0, 1.0);

    for (int i = 0; i < samples; ++i) {
        double t = ur(rng);
        // Mix log-like distribution near horizon and linear far out
        double rmin = r_plus + 1e-3 * std::max(1.0, M);
        double rmax = 50.0 * std::max(1.0, M);
        // Use a simple power interpolation to bias near rmin
        double r = rmin * std::pow(rmax / rmin, t);

        double th = (0.05 + 0.90 * ur(rng)) * BH_PI;
        // Steps: scale with magnitude to reduce subtractive cancellation
        double hr = std::max(1e-4 * std::max(1.0, r), 1e-6);
        double hth = std::max(1e-4, 1e-6);
        // Clamp theta steps to stay within (0, pi)
        if (th + 2*hth >= BH_PI) th = BH_PI - 2.5*hth;
        if (th - 2*hth <= 0.0)   th = 2.5*hth;

        // Analytic
        float gtt_r_a, gtphi_r_a, gphiphi_r_a, grr_r_a, gthth_r_a;
        float gtt_th_a, gtphi_th_a, gphiphi_th_a, grr_th_a, gthth_th_a;
        kerr_blocks_contra_derivs((float)M, (float)a, (float)r, (float)th,
                                  gtt_r_a, gtphi_r_a, gphiphi_r_a, grr_r_a, gthth_r_a,
                                  gtt_th_a, gtphi_th_a, gphiphi_th_a, grr_th_a, gthth_th_a);

        // Finite difference
        double gtt_r_f, gtphi_r_f, gphiphi_r_f, grr_r_f, gthth_r_f;
        double gtt_th_f, gtphi_th_f, gphiphi_th_f, grr_th_f, gthth_th_f;
        fd_derivs(M, a, r, th, hr, hth,
                  gtt_r_f, gtphi_r_f, gphiphi_r_f, grr_r_f, gthth_r_f,
                  gtt_th_f, gtphi_th_f, gphiphi_th_f, grr_th_f, gthth_th_f);

        // Track max rel errors
        st.max_rel_gtt_r     = std::max(st.max_rel_gtt_r,     rel_err(gtt_r_a,     gtt_r_f));
        st.max_rel_gtphi_r   = std::max(st.max_rel_gtphi_r,   rel_err(gtphi_r_a,   gtphi_r_f));
        st.max_rel_gphiphi_r = std::max(st.max_rel_gphiphi_r, rel_err(gphiphi_r_a, gphiphi_r_f));
        st.max_rel_grr_r     = std::max(st.max_rel_grr_r,     rel_err(grr_r_a,     grr_r_f));
        st.max_rel_gthth_r   = std::max(st.max_rel_gthth_r,   rel_err(gthth_r_a,   gthth_r_f));

        st.max_rel_gtt_th     = std::max(st.max_rel_gtt_th,     rel_err(gtt_th_a,     gtt_th_f));
        st.max_rel_gtphi_th   = std::max(st.max_rel_gtphi_th,   rel_err(gtphi_th_a,   gtphi_th_f));
        st.max_rel_gphiphi_th = std::max(st.max_rel_gphiphi_th, rel_err(gphiphi_th_a, gphiphi_th_f));
        st.max_rel_grr_th     = std::max(st.max_rel_grr_th,     rel_err(grr_th_a,     grr_th_f));
        st.max_rel_gthth_th   = std::max(st.max_rel_gthth_th,   rel_err(gthth_th_a,   gthth_th_f));

        // Count failures with mixed rel/abs tolerance
        auto ok = [&](double a1, double a2){
            double ra = std::abs(a1), rb = std::abs(a2);
            if (std::max(ra, rb) < 1e-7) {
                return std::abs(a1 - a2) <= abs_tol; // near-zero derivatives: use absolute tol
            } else {
                return rel_err(a1, a2) <= tol;
            }
        };
        int local_fails = 0;
        local_fails += !ok(gtt_r_a,     gtt_r_f);
        local_fails += !ok(gtphi_r_a,   gtphi_r_f);
        local_fails += !ok(gphiphi_r_a, gphiphi_r_f);
        local_fails += !ok(grr_r_a,     grr_r_f);
        local_fails += !ok(gthth_r_a,   gthth_r_f);
        local_fails += !ok(gtt_th_a,     gtt_th_f);
        local_fails += !ok(gtphi_th_a,   gtphi_th_f);
        local_fails += !ok(gphiphi_th_a, gphiphi_th_f);
        local_fails += !ok(grr_th_a,     grr_th_f);
        local_fails += !ok(gthth_th_a,   gthth_th_f);
        st.fails += (local_fails > 0);
    }
    return st;
}

int main(int argc, char** argv) {
    // Config via env or defaults
    int samples = 1000;
    double tol = 5e-3;     // float-precision target (adjust via BH_DERIV_TOL)
    double abs_tol = 5e-4; // absolute tol for near-zero derivatives
    unsigned seed = 12345u;
    double M = 1.0;      // use M=1 geometric units

    if (const char* s = std::getenv("BH_SAMPLES"))      samples = std::max(10, std::atoi(s));
    if (const char* s = std::getenv("BH_DERIV_TOL"))    tol = std::max(1e-9, std::atof(s));
    if (const char* s = std::getenv("BH_DERIV_ABS_TOL")) abs_tol = std::max(1e-12, std::atof(s));
    if (const char* s = std::getenv("BH_SEED"))         seed = (unsigned)std::strtoul(s, nullptr, 10);

    std::printf("[deriv-test] samples=%d tol=%.3e abs_tol=%.3e seed=%u M=%.3f\n", samples, tol, abs_tol, seed, M);
    std::vector<double> spins = {0.0, 0.5, 0.9};
    bool all_ok = true;
    for (double a_over_M : spins) {
        double a = a_over_M * M;
        RelErrStats st = run_suite(M, a, samples, tol, abs_tol, seed ^ (unsigned)std::llround(1000*a_over_M));
        std::printf("[spin a/M=%.2f] max |rel err| (r-derivs): gtt=%.2e gtphi=%.2e gphiphi=%.2e grr=%.2e gthth=%.2e\n",
                    a_over_M, st.max_rel_gtt_r, st.max_rel_gtphi_r, st.max_rel_gphiphi_r, st.max_rel_grr_r, st.max_rel_gthth_r);
        std::printf("[spin a/M=%.2f] max |rel err| (th-derivs): gtt=%.2e gtphi=%.2e gphiphi=%.2e grr=%.2e gthth=%.2e\n",
                    a_over_M, st.max_rel_gtt_th, st.max_rel_gtphi_th, st.max_rel_gphiphi_th, st.max_rel_grr_th, st.max_rel_gthth_th);
        std::printf("[spin a/M=%.2f] samples exceeding tol: %d / %d\n", a_over_M, st.fails, samples);
        all_ok = all_ok && (st.fails == 0);
    }

    if (all_ok) {
        std::puts("[pass] Kerr analytic derivative self-test within tolerance.");
        return 0;
    } else {
        std::puts("[fail] Some samples exceeded derivative tolerance; inspect logs.");
        return 1;
    }
}
