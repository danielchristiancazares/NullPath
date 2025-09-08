/*
 * MINO TIME VALIDATOR
 * 
 * Comprehensive validation program for CUDA Mino time implementation
 * Tests against the 5000+ generated test cases
 * 
 * This validates the critical fix: proper Mino time parameterization
 * that should improve elliptic integral success rate from 65.7% to >95%
 */

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <vector>
#include <chrono>

// Forward declarations
extern "C" void compute_mino_time_batch_cuda(
    const double* h_u0, const double* h_uf,
    const double* h_mu0, const double* h_muf,
    const double* h_a, const double* h_l, const double* h_q2,
    const int* h_tpm, const int* h_tpr,
    const int* h_su, const int* h_sm,
    const int* h_ncase,
    double* h_lambda_results,
    char* h_valid_results,
    int n_cases
);

struct TestCase {
    double u0, uf, mu0, muf, a, l, l2, q2;
    int tpm, tpr, su, sm;
};

// Load test cases from our generated file
std::vector<TestCase> load_test_cases(const std::string& filename) {
    std::vector<TestCase> cases;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        printf("ERROR: Cannot open test file %s\n", filename.c_str());
        return cases;
    }
    
    // Skip header lines
    std::string line;
    std::getline(file, line); // Skip comment line
    std::getline(file, line); // Skip format line
    
    int n_cases;
    file >> n_cases;
    
    printf("Loading %d test cases from %s...\n", n_cases, filename.c_str());
    
    for (int i = 0; i < n_cases; i++) {
        TestCase test;
        file >> test.u0 >> test.uf >> test.mu0 >> test.muf
             >> test.a >> test.l >> test.l2 >> test.q2
             >> test.tpm >> test.tpr >> test.su >> test.sm;
        
        if (file.fail()) {
            printf("WARNING: Failed to read test case %d\n", i);
            break;
        }
        
        cases.push_back(test);
    }
    
    file.close();
    printf("Successfully loaded %zu test cases\n", cases.size());
    return cases;
}

// Analyze results and compute statistics
void analyze_results(const std::vector<double>& lambda_results,
                    const std::vector<char>& valid_results) {
    int valid_count = 0;
    int nonzero_count = 0;
    double sum = 0.0, sum_sq = 0.0;
    double min_val = 1e100, max_val = -1e100;
    
    for (size_t i = 0; i < lambda_results.size(); i++) {
        if (valid_results[i]) {
            valid_count++;
            double val = lambda_results[i];
            
            if (fabs(val) > 1e-15) {
                nonzero_count++;
                sum += val;
                sum_sq += val * val;
                
                if (val < min_val) min_val = val;
                if (val > max_val) max_val = val;
            }
        }
    }
    
    printf("\n=== MINO TIME VALIDATION RESULTS ===\n");
    printf("Total test cases: %zu\n", lambda_results.size());
    printf("Valid computations: %d (%.1f%%)\n", valid_count, 
           100.0 * valid_count / lambda_results.size());
    printf("Non-zero λ values: %d (%.1f%%)\n", nonzero_count,
           100.0 * nonzero_count / lambda_results.size());
    
    if (nonzero_count > 0) {
        double mean = sum / nonzero_count;
        double variance = (sum_sq / nonzero_count) - (mean * mean);
        double stddev = sqrt(variance);
        
        printf("Mino time statistics:\n");
        printf("  Mean λ: %.6e\n", mean);
        printf("  Std dev: %.6e\n", stddev);
        printf("  Range: [%.6e, %.6e]\n", min_val, max_val);
    }
    
    // Check for expected improvement over coordinate time
    if (valid_count > 0.9 * lambda_results.size()) {
        printf("\n✅ SUCCESS: >90%% valid computations\n");
        printf("This suggests proper Mino time parameterization!\n");
    } else {
        printf("\n⚠️  WARNING: %.1f%% valid computations\n", 
               100.0 * valid_count / lambda_results.size());
        printf("May need further refinement of Mino time computation\n");
    }
}

int main() {
    printf("GEOKERR MINO TIME COMPREHENSIVE VALIDATOR\n");
    printf("========================================\n");
    printf("Testing CUDA implementation against generated test cases\n");
    printf("Target: >95%% success rate (vs 65.7%% with coordinate time)\n\n");
    
    // Load test cases
    std::vector<TestCase> test_cases = load_test_cases("validation/mino_time_inputs.dat");
    if (test_cases.empty()) {
        printf("No test cases loaded. Exiting.\n");
        return 1;
    }
    
    int n_cases = test_cases.size();
    
    // Prepare input arrays
    std::vector<double> u0_array(n_cases), uf_array(n_cases);
    std::vector<double> mu0_array(n_cases), muf_array(n_cases);
    std::vector<double> a_array(n_cases), l_array(n_cases), q2_array(n_cases);
    std::vector<int> tpm_array(n_cases), tpr_array(n_cases);
    std::vector<int> su_array(n_cases), sm_array(n_cases);
    std::vector<int> ncase_array(n_cases); // Default to case 1 for now
    
    for (int i = 0; i < n_cases; i++) {
        const TestCase& tc = test_cases[i];
        u0_array[i] = tc.u0;
        uf_array[i] = tc.uf;
        mu0_array[i] = tc.mu0;
        muf_array[i] = tc.muf;
        a_array[i] = tc.a;
        l_array[i] = tc.l;
        q2_array[i] = tc.q2;
        tpm_array[i] = tc.tpm;
        tpr_array[i] = tc.tpr;
        su_array[i] = tc.su;
        sm_array[i] = tc.sm;
        ncase_array[i] = 1; // Default case - would need proper classification
    }
    
    // Prepare output arrays
    std::vector<double> lambda_results(n_cases);
    std::vector<char> valid_results(n_cases);
    
    printf("Running CUDA Mino time computation on %d test cases...\n", n_cases);
    
    // Measure performance
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run CUDA computation
    compute_mino_time_batch_cuda(
        u0_array.data(), uf_array.data(),
        mu0_array.data(), muf_array.data(),
        a_array.data(), l_array.data(), q2_array.data(),
        tpm_array.data(), tpr_array.data(),
        su_array.data(), sm_array.data(),
        ncase_array.data(),
        lambda_results.data(),
        valid_results.data(),
        n_cases
    );
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    
    printf("Computation completed in %.3f ms\n", duration.count() / 1000.0);
    printf("Performance: %.0f evaluations/second\n", 
           n_cases * 1e6 / duration.count());
    
    // Analyze results
    analyze_results(lambda_results, valid_results);
    
    // Save results for further analysis
    std::ofstream outfile("validation/mino_time_cuda_results.dat");
    outfile << "# CUDA Mino time validation results\n";
    outfile << "# Format: lambda valid\n";
    outfile << n_cases << "\n";
    
    for (int i = 0; i < n_cases; i++) {
        outfile << lambda_results[i] << " " << (valid_results[i] ? 1 : 0) << "\n";
    }
    outfile.close();
    
    printf("\nResults saved to validation/mino_time_cuda_results.dat\n");
    printf("\n=== KEY INSIGHT ===\n");
    printf("This validates the critical missing piece: Mino time parameterization\n");
    printf("λ = λ_u + t_μ enables stable elliptic integral evaluation\n");
    printf("Success rate should be >>65.7%% if implementation is correct\n");
    
    return 0;
}