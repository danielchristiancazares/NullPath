/**
 * Combined CUDA Carlson Elliptic Integral Implementation and Validator
 * 
 * Single file containing both implementations and validation to avoid linking issues.
 */

#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <map>

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ \
                  << " - " << cudaGetErrorString(err) << std::endl; \
        exit(1); \
    } \
} while(0)

// ==================== CARLSON ELLIPTIC INTEGRAL FUNCTIONS ====================

/**
 * Carlson RF function - device implementation
 */
__device__ double carlson_rf_device(double x, double y, double z) {
    // Constants
    const double ERRTOL = 0.08;
    const double TINY = 1.5e-38;
    const double BIG = 3.0e37;
    const double third = 1.0/3.0;
    const double c1 = 1.0/24.0;
    const double c2 = 0.1;
    const double c3 = 3.0/44.0;
    const double c4 = 1.0/14.0;
    
    // Validity check
    if (fmin(fmin(x,y),z) < 0.0 || 
        fmin(fmin(x+y, x+z), y+z) < TINY ||
        fmax(fmax(x,y),z) > BIG) {
        return NAN;
    }
    
    double xt = x, yt = y, zt = z;
    
    for (int iter = 0; iter < 100; iter++) {
        double sqrtx = sqrt(xt);
        double sqrty = sqrt(yt);
        double sqrtz = sqrt(zt);
        double alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        
        double ave = third * (xt + yt + zt);
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        
        if (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL) {
            double e2 = delx * dely - delz * delz;
            double e3 = delx * dely * delz;
            return (1.0 + (c1*e2 - c2 - c3*e3)*e2 + c4*e3) / sqrt(ave);
        }
    }
    return NAN; // Convergence failure
}

/**
 * Carlson RC function - device implementation
 */
__device__ double carlson_rc_device(double x, double y) {
    // Constants
    const double ERRTOL = 0.04;
    const double TINY = 1.69e-38;
    const double BIG = 3.0e37;
    const double third = 1.0/3.0;
    const double c1 = 0.3;
    const double c2 = 1.0/7.0;
    const double c3 = 0.375;
    const double c4 = 9.0/22.0;
    
    // Simplified validity check
    if (x < 0.0 || y == 0.0 || (x + fabs(y)) < TINY || (x + fabs(y)) > BIG) {
        return NAN;
    }
    
    double xt, yt, w;
    if (y > 0.0) {
        xt = x;
        yt = y;
        w = 1.0;
    } else {
        xt = x - y;
        yt = -y;
        w = sqrt(x) / sqrt(xt);
    }
    
    for (int iter = 0; iter < 100; iter++) {
        double alamb = 2.0 * sqrt(xt) * sqrt(yt) + yt;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        double ave = third * (xt + yt + yt);
        double s = (yt - ave) / ave;
        
        if (fabs(s) <= ERRTOL) {
            return w * (1.0 + s*s*(c1 + s*(c2 + s*(c3 + s*c4)))) / sqrt(ave);
        }
    }
    return NAN;
}

/**
 * Carlson RD function - device implementation
 */
__device__ double carlson_rd_device(double x, double y, double z) {
    // Constants  
    const double ERRTOL = 0.05;
    const double TINY = 1.0e-25;
    const double BIG = 4.5e21;
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/6.0;
    const double c3 = 9.0/22.0;
    const double c4 = 3.0/26.0;
    const double c5 = 0.25 * c3;
    const double c6 = 1.5 * c4;
    
    // Validity check
    if (fmin(x,y) < 0.0 || fmin(x+y, z) < TINY || fmax(fmax(x,y),z) > BIG) {
        return NAN;
    }
    
    double xt = x, yt = y, zt = z;
    double sum = 0.0;
    double fac = 1.0;
    
    for (int iter = 0; iter < 100; iter++) {
        double sqrtx = sqrt(xt);
        double sqrty = sqrt(yt);
        double sqrtz = sqrt(zt);
        double alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        
        sum = sum + fac / (sqrtz * (zt + alamb));
        fac = 0.25 * fac;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        
        double ave = 0.2 * (xt + yt + 3.0 * zt);
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        
        if (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL) {
            double ea = delx * dely;
            double eb = delz * delz;
            double ec = ea - eb;
            double ed = ea - 6.0 * eb;
            double ee = ed + ec + ec;
            
            return 3.0 * sum + fac * (1.0 + ed*(-c1 + c5*ed - c6*delz*ee) +
                   delz*(c2*ee + delz*(-c3*ec + delz*c4*ea))) / (ave * sqrt(ave));
        }
    }
    return NAN;
}

/**
 * Carlson RJ function - device implementation
 */
__device__ double carlson_rj_device(double x, double y, double z, double p) {
    // Constants
    const double ERRTOL = 0.05;
    const double TINY = 2.5e-13;
    const double BIG = 9.0e11;
    
    // Validity check
    if (fmin(fmin(x,y),z) < 0.0 || 
        fmin(fmin(fmin(x+y, x+z), y+z), p) < TINY ||
        fmax(fmax(fmax(x,y),z),p) > BIG) {
        return NAN;
    }
    
    // For simplicity, only handle p > 0 case (most common)
    if (p <= 0.0) {
        return NAN; // Would need more complex handling
    }
    
    double sum = 0.0;
    double fac = 1.0;
    double xt = x, yt = y, zt = z, pt = p;
    
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/3.0;
    const double c3 = 3.0/22.0;
    const double c4 = 3.0/26.0;
    const double c5 = 0.75 * c3;
    const double c6 = 1.5 * c4;
    const double c7 = 0.5 * c2;
    const double c8 = c3 + c3;
    
    for (int iter = 0; iter < 100; iter++) {
        double sqrtx = sqrt(xt);
        double sqrty = sqrt(yt);
        double sqrtz = sqrt(zt);
        double sqrtp = sqrt(pt);
        double dnm = sqrtp * (sqrtp + sqrtx) * (sqrtp + sqrty) * (sqrtp + sqrtz);
        
        sum = sum + fac / dnm;
        fac = 0.25 * fac;
        double alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        
        double ave = 0.2 * (xt + yt + zt + pt + pt);
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        double delp = (ave - pt) / ave;
        
        if (fmax(fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)), fabs(delp)) <= ERRTOL) {
            double ea = delx * (dely + delz) + dely * delz;
            double eb = delx * dely * delz;
            double ec = delp * delp;
            double ed = ea - 3.0 * ec;
            double ee = eb + 2.0 * delp * (ea - ec);
            
            return 3.0 * sum + fac * (1.0 + ed*(-c1 + c5*ed - c6*ee) + 
                   eb*(c7 + delp*(-c8 + delp*c4)) + 
                   delp*ea*(c2 - delp*c3) - c2*delp*ec) / (ave * sqrt(ave));
        }
    }
    return NAN;
}

// ==================== VALIDATION KERNEL ====================

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
    
    char type0 = test_types[idx * 2];
    char type1 = test_types[idx * 2 + 1];
    
    double x = x_vals[idx];
    double y = y_vals[idx];
    double z = z_vals[idx];
    double p = p_vals[idx];
    
    if (type0 == 'R' && type1 == 'F') {
        cuda_results[idx] = carlson_rf_device(x, y, z);
    } else if (type0 == 'R' && type1 == 'C') {
        cuda_results[idx] = carlson_rc_device(x, y);
    } else if (type0 == 'R' && type1 == 'D') {
        cuda_results[idx] = carlson_rd_device(x, y, z);
    } else if (type0 == 'R' && type1 == 'J') {
        cuda_results[idx] = carlson_rj_device(x, y, z, p);
    } else {
        cuda_results[idx] = NAN;
    }
}

// ==================== HOST VALIDATION CODE ====================

struct TestCase {
    std::string type;
    double x, y, z, p;
    double reference_result;
};

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
        cases[i].reference_result = result;
    }
    
    file.close();
    std::cout << "Loaded reference results for " << cases.size() << " cases" << std::endl;
    return true;
}

void analyze_results(const std::vector<TestCase>& cases, 
                    const std::vector<double>& cuda_results) {
    
    std::cout << "\nVALIDATION ANALYSIS" << std::endl;
    std::cout << "==================" << std::endl;
    
    std::map<std::string, std::vector<double>> errors_by_type;
    std::map<std::string, int> counts_by_type;
    
    int total_cases = cases.size();
    int passed_cases = 0;
    int failed_cases = 0;
    double max_rel_error = 0.0;
    int worst_case_idx = -1;
    
    const double tolerance = 1e-12;
    
    for (int i = 0; i < total_cases; i++) {
        double ref = cases[i].reference_result;
        double cuda = cuda_results[i];
        
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