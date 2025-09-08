/**
 * CUDA Carlson Elliptic Integral Validator
 * 
 * Loads test cases and reference data, runs CUDA implementations,
 * and validates results against FORTRAN reference with detailed statistics.
 */

#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <map>

// Forward declarations - functions are defined in carlson_elliptic.cu
__device__ double carlson_rf(double x, double y, double z);
__device__ double carlson_rc(double x, double y);
__device__ double carlson_rd(double x, double y, double z);
__device__ double carlson_rj(double x, double y, double z, double p);

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ \
                  << " - " << cudaGetErrorString(err) << std::endl; \
        exit(1); \
    } \
} while(0)

struct TestCase {
    std::string type;
    double x, y, z, p;
    double reference_result;
};

/**
 * CUDA kernel to batch process elliptic integral validation
 */
__global__ void validate_elliptic_kernel(
    const char* test_types, 
    const double* x_vals, const double* y_vals, 
    const double* z_vals, const double* p_vals,
    double* cuda_results, int n_tests
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_tests) return;
    
    // Extract test type (2 characters)
    char type0 = test_types[idx * 2];
    char type1 = test_types[idx * 2 + 1];
    
    double x = x_vals[idx];
    double y = y_vals[idx];
    double z = z_vals[idx];
    double p = p_vals[idx];
    
    if (type0 == 'R' && type1 == 'F') {
        cuda_results[idx] = carlson_rf(x, y, z);
    } else if (type0 == 'R' && type1 == 'C') {
        cuda_results[idx] = carlson_rc(x, y);
    } else if (type0 == 'R' && type1 == 'D') {
        cuda_results[idx] = carlson_rd(x, y, z);
    } else if (type0 == 'R' && type1 == 'J') {
        cuda_results[idx] = carlson_rj(x, y, z, p);
    } else {
        cuda_results[idx] = NAN; // Unknown test type
    }
}

/**
 * Load test cases from FORTRAN input file
 */
std::vector<TestCase> load_input_cases(const std::string& filename) {
    std::vector<TestCase> cases;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open input file " << filename << std::endl;
        return cases;
    }
    
    int n_cases;
    file >> n_cases;
    
    std::cout << "Loading " << n_cases << " test cases from " << filename << std::endl;
    
    for (int i = 0; i < n_cases; i++) {
        TestCase tc;
        file >> tc.type >> tc.x >> tc.y >> tc.z >> tc.p;
        cases.push_back(tc);
    }
    
    file.close();
    std::cout << "Loaded " << cases.size() << " test cases" << std::endl;
    return cases;
}

/**
 * Load reference results from FORTRAN output file
 */
bool load_reference_results(const std::string& filename, std::vector<TestCase>& cases) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open reference file " << filename << std::endl;
        return false;
    }
    
    int n_cases;
    file >> n_cases;
    
    if (n_cases != cases.size()) {
        std::cerr << "Error: Reference file has " << n_cases 
                  << " cases but input has " << cases.size() << std::endl;
        return false;
    }
    
    std::cout << "Loading reference results from " << filename << std::endl;
    
    for (int i = 0; i < n_cases; i++) {
        std::string type;
        double x, y, z, p, result;
        file >> type >> x >> y >> z >> p >> result;
        
        // Verify this matches the input case
        if (type != cases[i].type || 
            fabs(x - cases[i].x) > 1e-15 ||
            fabs(y - cases[i].y) > 1e-15 ||
            fabs(z - cases[i].z) > 1e-15 ||
            fabs(p - cases[i].p) > 1e-15) {
            std::cerr << "Error: Case " << i << " mismatch between input and reference" << std::endl;
            return false;
        }
        
        cases[i].reference_result = result;
    }
    
    file.close();
    std::cout << "Loaded reference results for " << cases.size() << " cases" << std::endl;
    return true;
}

/**
 * Analyze validation results and print statistics
 */
void analyze_results(const std::vector<TestCase>& cases, 
                    const std::vector<double>& cuda_results) {
    
    std::cout << "\nVALIDATION ANALYSIS" << std::endl;
    std::cout << "==================" << std::endl;
    
    // Statistics by function type
    std::map<std::string, std::vector<double>> errors_by_type;
    std::map<std::string, int> counts_by_type;
    
    int total_cases = cases.size();
    int passed_cases = 0;
    int failed_cases = 0;
    double max_rel_error = 0.0;
    int worst_case_idx = -1;
    
    const double tolerance = 1e-12; // Target accuracy
    
    for (int i = 0; i < total_cases; i++) {
        double ref = cases[i].reference_result;
        double cuda = cuda_results[i];
        
        // Skip NaN results (likely invalid inputs)
        if (std::isnan(ref) || std::isnan(cuda)) {
            continue;
        }
        
        double abs_error = fabs(cuda - ref);
        double rel_error = (fabs(ref) > 1e-16) ? abs_error / fabs(ref) : abs_error;
        
        errors_by_type[cases[i].type].push_back(rel_error);
        counts_by_type[cases[i].type]++;
        
        if (rel_error > max_rel_error) {
            max_rel_error = rel_error;
            worst_case_idx = i;
        }
        
        if (rel_error <= tolerance) {
            passed_cases++;
        } else {
            failed_cases++;
            
            // Print details for significant failures
            if (rel_error > 1e-10) {
                std::cout << "FAIL [" << i << "]: " << cases[i].type 
                          << "(" << cases[i].x << ", " << cases[i].y 
                          << ", " << cases[i].z << ", " << cases[i].p << ")" << std::endl;
                std::cout << "  Reference: " << std::scientific << std::setprecision(12) << ref << std::endl;
                std::cout << "  CUDA:      " << std::scientific << std::setprecision(12) << cuda << std::endl;
                std::cout << "  Rel Error: " << std::scientific << std::setprecision(6) << rel_error << std::endl;
            }
        }
    }
    
    // Overall statistics
    std::cout << "\nOVERALL RESULTS:" << std::endl;
    std::cout << "  Total cases:     " << total_cases << std::endl;
    std::cout << "  Passed:          " << passed_cases << " (" 
              << std::fixed << std::setprecision(2) 
              << 100.0 * passed_cases / total_cases << "%)" << std::endl;
    std::cout << "  Failed:          " << failed_cases << " (" 
              << std::fixed << std::setprecision(2)
              << 100.0 * failed_cases / total_cases << "%)" << std::endl;
    std::cout << "  Max rel error:   " << std::scientific << std::setprecision(6) 
              << max_rel_error << std::endl;
    
    if (worst_case_idx >= 0) {
        std::cout << "\nWORST CASE [" << worst_case_idx << "]:" << std::endl;
        std::cout << "  Function: " << cases[worst_case_idx].type << std::endl;
        std::cout << "  Args: (" << cases[worst_case_idx].x << ", " 
                  << cases[worst_case_idx].y << ", " << cases[worst_case_idx].z 
                  << ", " << cases[worst_case_idx].p << ")" << std::endl;
        std::cout << "  Reference: " << std::setprecision(12) 
                  << cases[worst_case_idx].reference_result << std::endl;
        std::cout << "  CUDA:      " << std::setprecision(12) 
                  << cuda_results[worst_case_idx] << std::endl;
    }
    
    // Per-function statistics
    std::cout << "\nPER-FUNCTION ANALYSIS:" << std::endl;
    for (const auto& pair : errors_by_type) {
        const std::string& func = pair.first;
        const std::vector<double>& errors = pair.second;
        
        if (errors.empty()) continue;
        
        double mean_error = 0.0;
        double max_error = 0.0;
        for (double err : errors) {
            mean_error += err;
            max_error = std::max(max_error, err);
        }
        mean_error /= errors.size();
        
        std::cout << "  " << func << ": " << counts_by_type[func] << " cases, "
                  << "mean=" << std::scientific << std::setprecision(3) << mean_error
                  << ", max=" << std::scientific << std::setprecision(3) << max_error
                  << std::endl;
    }
}

int main() {
    std::cout << "CUDA CARLSON ELLIPTIC INTEGRAL VALIDATOR" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Load test cases
    std::vector<TestCase> cases = load_input_cases("../validation/test_data/elliptic_inputs.dat");
    if (cases.empty()) {
        std::cerr << "Error: No test cases loaded" << std::endl;
        return 1;
    }
    
    // Load reference results
    if (!load_reference_results("../validation/test_data/elliptic_reference.dat", cases)) {
        std::cerr << "Error: Failed to load reference results" << std::endl;
        return 1;
    }
    
    int n_cases = cases.size();
    
    // Prepare GPU data
    std::vector<char> types_flat(n_cases * 2);
    std::vector<double> x_vals(n_cases), y_vals(n_cases), z_vals(n_cases), p_vals(n_cases);
    
    for (int i = 0; i < n_cases; i++) {
        types_flat[i*2] = cases[i].type[0];
        types_flat[i*2+1] = cases[i].type[1];
        x_vals[i] = cases[i].x;
        y_vals[i] = cases[i].y;
        z_vals[i] = cases[i].z;
        p_vals[i] = cases[i].p;
    }
    
    // Allocate GPU memory
    char *d_types;
    double *d_x, *d_y, *d_z, *d_p, *d_results;
    
    CUDA_CHECK(cudaMalloc(&d_types, n_cases * 2));
    CUDA_CHECK(cudaMalloc(&d_x, n_cases * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_y, n_cases * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_z, n_cases * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_p, n_cases * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_results, n_cases * sizeof(double)));
    
    // Copy data to GPU
    CUDA_CHECK(cudaMemcpy(d_types, types_flat.data(), n_cases * 2, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_x, x_vals.data(), n_cases * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_y, y_vals.data(), n_cases * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z, z_vals.data(), n_cases * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_p, p_vals.data(), n_cases * sizeof(double), cudaMemcpyHostToDevice));
    
    // Launch kernel
    int block_size = 256;
    int grid_size = (n_cases + block_size - 1) / block_size;
    
    std::cout << "\nRunning CUDA validation on " << n_cases << " cases..." << std::endl;
    std::cout << "Grid size: " << grid_size << ", Block size: " << block_size << std::endl;
    
    // Time the kernel execution
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    
    CUDA_CHECK(cudaEventRecord(start));
    validate_elliptic_kernel<<<grid_size, block_size>>>(
        d_types, d_x, d_y, d_z, d_p, d_results, n_cases);
    CUDA_CHECK(cudaEventRecord(stop));
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    float gpu_time_ms;
    CUDA_CHECK(cudaEventElapsedTime(&gpu_time_ms, start, stop));
    
    std::cout << "GPU computation completed in " << gpu_time_ms << " ms" << std::endl;
    std::cout << "Throughput: " << std::fixed << std::setprecision(0) 
              << n_cases / (gpu_time_ms / 1000.0) << " evaluations/second" << std::endl;
    
    // Copy results back
    std::vector<double> cuda_results(n_cases);
    CUDA_CHECK(cudaMemcpy(cuda_results.data(), d_results, n_cases * sizeof(double), cudaMemcpyDeviceToHost));
    
    // Analyze and report results
    analyze_results(cases, cuda_results);
    
    // Cleanup
    cudaFree(d_types);
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(d_p);
    cudaFree(d_results);
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    return 0;
}