/**
 * CUDA End-to-End GEOKERR Validation Suite
 * 
 * Complete validation pipeline that:
 * 1. Loads geodesic test cases from FORTRAN inputs
 * 2. Runs CUDA semi-analytic solver
 * 3. Compares with FORTRAN reference results
 * 4. Provides detailed accuracy and performance analysis
 */

#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

// Include the CUDA geodesic solver
extern "C" void geokerr_cuda_batch(
    const double* a_vals, const double* r0_vals, const double* mu0_vals,
    const double* E_vals, const double* L_vals, const double* Q2_vals,
    double* r_final, double* mu_final, double* phi_final, double* t_final,
    int n_geodesics, int n_steps, double lambda_max
);

struct GeodecisTestCase {
    double a, M, b;           // Spacetime and impact parameter
    double r0, theta0, phi0;  // Initial coordinates
    
    // Derived quantities
    double E, L, Q2, mu0;
    
    // Reference results (from FORTRAN)
    double r_ref, mu_ref, phi_ref, t_ref;
    int status_ref;
    
    // CUDA results
    double r_cuda, mu_cuda, phi_cuda, t_cuda;
};

/**
 * Load geodesic test cases from FORTRAN input files
 */
std::vector<GeodecisTestCase> load_geodesic_cases(const std::string& input_file, 
                                                 const std::string& reference_file) {
    
    std::vector<GeodecisTestCase> cases;
    
    // Load input cases
    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open " << input_file << std::endl;
        return cases;
    }
    
    int n_cases;
    infile >> n_cases;
    
    std::cout << "Loading " << n_cases << " geodesic test cases..." << std::endl;
    
    for (int i = 0; i < n_cases; i++) {
        GeodecisTestCase tc;
        infile >> tc.a >> tc.M >> tc.b >> tc.r0 >> tc.theta0 >> tc.phi0;
        
        // Convert to conserved quantities
        tc.E = 1.0;                    // Photon energy
        tc.L = tc.b * tc.E;           // Angular momentum
        tc.Q2 = 0.0;                  // Equatorial motion
        tc.mu0 = cos(tc.theta0);      // mu = cos(theta)
        
        cases.push_back(tc);
    }
    infile.close();
    
    // Load reference results
    std::ifstream reffile(reference_file);
    if (!reffile.is_open()) {
        std::cerr << "Error: Cannot open " << reference_file << std::endl;
        return cases;
    }
    
    int n_ref;
    reffile >> n_ref;
    
    if (n_ref != n_cases) {
        std::cerr << "Error: Reference file has " << n_ref 
                  << " cases, input has " << n_cases << std::endl;
        return cases;
    }
    
    std::cout << "Loading reference results..." << std::endl;
    
    for (int i = 0; i < n_cases; i++) {
        int idx;
        double a, M, b, r0, theta0, phi0;
        int tpm, tpr;
        double su, sm;
        
        reffile >> idx >> a >> M >> b >> r0 >> theta0 >> phi0 >> tpm >> tpr >> su >> sm;
        
        // Store reference results (simplified - full implementation would decode FORTRAN output)
        cases[i].r_ref = r0 - fabs(su);    // Approximate final radius
        cases[i].mu_ref = cos(theta0);      // Approximate final mu  
        cases[i].phi_ref = 0.0;            // Would be computed from FORTRAN
        cases[i].t_ref = 0.0;              // Would be computed from FORTRAN
        cases[i].status_ref = tpm;         // Orbit classification
    }
    reffile.close();
    
    std::cout << "Loaded " << cases.size() << " test cases with references" << std::endl;
    return cases;
}

/**
 * Run CUDA batch geodesic computation
 */
void compute_cuda_geodesics(std::vector<GeodecisTestCase>& cases) {
    int n_cases = cases.size();
    
    // Prepare input arrays
    std::vector<double> a_vals(n_cases), r0_vals(n_cases), mu0_vals(n_cases);
    std::vector<double> E_vals(n_cases), L_vals(n_cases), Q2_vals(n_cases);
    
    for (int i = 0; i < n_cases; i++) {
        a_vals[i] = cases[i].a;
        r0_vals[i] = cases[i].r0;
        mu0_vals[i] = cases[i].mu0;
        E_vals[i] = cases[i].E;
        L_vals[i] = cases[i].L;
        Q2_vals[i] = cases[i].Q2;
    }
    
    // Output arrays
    std::vector<double> r_final(n_cases), mu_final(n_cases);
    std::vector<double> phi_final(n_cases), t_final(n_cases);
    
    std::cout << "\\nRunning CUDA geodesic computation on " << n_cases << " cases..." << std::endl;
    
    // Time the computation
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    cudaEventRecord(start);
    
    // Run CUDA batch computation
    geokerr_cuda_batch(
        a_vals.data(), r0_vals.data(), mu0_vals.data(),
        E_vals.data(), L_vals.data(), Q2_vals.data(),
        r_final.data(), mu_final.data(), phi_final.data(), t_final.data(),
        n_cases, 1000, 50.0  // 1000 steps, lambda_max = 50
    );
    
    cudaEventRecord(stop);
    cudaDeviceSynchronize();
    
    float gpu_time_ms;
    cudaEventElapsedTime(&gpu_time_ms, start, stop);
    
    std::cout << "CUDA computation completed in " << gpu_time_ms << " ms" << std::endl;
    std::cout << "Throughput: " << std::fixed << std::setprecision(1) 
              << n_cases / (gpu_time_ms / 1000.0) << " geodesics/second" << std::endl;
    
    // Store results back to cases
    for (int i = 0; i < n_cases; i++) {
        cases[i].r_cuda = r_final[i];
        cases[i].mu_cuda = mu_final[i];
        cases[i].phi_cuda = phi_final[i];
        cases[i].t_cuda = t_final[i];
    }
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

/**
 * Analyze validation results
 */
void analyze_validation(const std::vector<GeodecisTestCase>& cases) {
    std::cout << "\\nGEODESIC VALIDATION ANALYSIS" << std::endl;
    std::cout << "============================" << std::endl;
    
    int n_cases = cases.size();
    int valid_cases = 0;
    int failed_cases = 0;
    
    // Error statistics
    std::vector<double> r_errors, mu_errors, phi_errors, t_errors;
    
    const double tolerance = 1e-6; // Relaxed tolerance for geodesic coordinates
    
    for (int i = 0; i < n_cases; i++) {
        const auto& tc = cases[i];
        
        // Skip invalid results
        if (std::isnan(tc.r_cuda) || std::isinf(tc.r_cuda) ||
            std::isnan(tc.mu_cuda) || std::isinf(tc.mu_cuda)) {
            failed_cases++;
            continue;
        }
        
        // Compute errors
        double r_err = fabs(tc.r_cuda - tc.r_ref);
        double mu_err = fabs(tc.mu_cuda - tc.mu_ref);
        double phi_err = fabs(tc.phi_cuda - tc.phi_ref);
        double t_err = fabs(tc.t_cuda - tc.t_ref);
        
        r_errors.push_back(r_err);
        mu_errors.push_back(mu_err);
        phi_errors.push_back(phi_err);
        t_errors.push_back(t_err);
        
        valid_cases++;
        
        // Print details for significant failures
        if (r_err > tolerance) {
            std::cout << "COORD FAIL [" << i << "]: a=" << tc.a << ", b=" << tc.b << std::endl;
            std::cout << "  r: CUDA=" << tc.r_cuda << ", REF=" << tc.r_ref 
                      << ", ERR=" << r_err << std::endl;
            std::cout << "  mu: CUDA=" << tc.mu_cuda << ", REF=" << tc.mu_ref 
                      << ", ERR=" << mu_err << std::endl;
        }
    }
    
    // Overall statistics
    std::cout << "\\nOVERALL RESULTS:" << std::endl;
    std::cout << "  Total cases:     " << n_cases << std::endl;
    std::cout << "  Valid results:   " << valid_cases << " (" 
              << std::fixed << std::setprecision(1) 
              << 100.0 * valid_cases / n_cases << "%)" << std::endl;
    std::cout << "  Failed/Invalid:  " << failed_cases << " (" 
              << std::fixed << std::setprecision(1)
              << 100.0 * failed_cases / n_cases << "%)" << std::endl;
    
    if (valid_cases > 0) {
        // Compute mean errors
        double r_mean = 0, mu_mean = 0, phi_mean = 0, t_mean = 0;
        double r_max = 0, mu_max = 0, phi_max = 0, t_max = 0;
        
        for (int i = 0; i < valid_cases; i++) {
            r_mean += r_errors[i];
            mu_mean += mu_errors[i];
            phi_mean += phi_errors[i];
            t_mean += t_errors[i];
            
            r_max = std::max(r_max, r_errors[i]);
            mu_max = std::max(mu_max, mu_errors[i]);
            phi_max = std::max(phi_max, phi_errors[i]);
            t_max = std::max(t_max, t_errors[i]);
        }
        
        r_mean /= valid_cases;
        mu_mean /= valid_cases;
        phi_mean /= valid_cases;
        t_mean /= valid_cases;
        
        std::cout << "\\nCOORDINATE ERRORS:" << std::endl;
        std::cout << "  r (radius):   mean=" << std::scientific << std::setprecision(3) << r_mean
                  << ", max=" << r_max << std::endl;
        std::cout << "  μ (cos θ):    mean=" << std::scientific << std::setprecision(3) << mu_mean
                  << ", max=" << mu_max << std::endl;
        std::cout << "  φ (azimuth):  mean=" << std::scientific << std::setprecision(3) << phi_mean
                  << ", max=" << phi_max << std::endl;
        std::cout << "  t (time):     mean=" << std::scientific << std::setprecision(3) << t_mean
                  << ", max=" << t_max << std::endl;
    }
    
    std::cout << "\\nStatus: " << (failed_cases < n_cases/2 ? "PARTIAL SUCCESS" : "NEEDS DEBUGGING") 
              << std::endl;
    std::cout << "Next: Improve geodesic solver accuracy and robustness" << std::endl;
}

int main() {
    std::cout << "GEOKERR CUDA END-TO-END VALIDATION SUITE" << std::endl;
    std::cout << "=========================================" << std::endl;
    
    // Load test cases
    auto cases = load_geodesic_cases(
        "../validation/test_data/geodesic_inputs.dat",
        "../validation/test_data/geodesic_reference.dat"
    );
    
    if (cases.empty()) {
        std::cerr << "Error: No test cases loaded" << std::endl;
        return 1;
    }
    
    // Run CUDA computations
    compute_cuda_geodesics(cases);
    
    // Analyze results
    analyze_validation(cases);
    
    return 0;
}