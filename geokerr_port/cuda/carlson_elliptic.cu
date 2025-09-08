/**
 * CUDA Implementation of Carlson Elliptic Integral Functions
 * 
 * Based on Carlson (1995): "Numerical computation of real or complex elliptic integrals"
 * Ported from FORTRAN geokerr implementation by Dexter & Agol (2009)
 * 
 * Functions:
 * - carlson_rf(x, y, z): Elliptic integral of the first kind
 * - carlson_rc(x, y):    Degenerate elliptic integral  
 * - carlson_rd(x, y, z): Elliptic integral of the second kind
 * - carlson_rj(x, y, z, p): Elliptic integral of the third kind
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <iostream>

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ \
                  << " - " << cudaGetErrorString(err) << std::endl; \
        exit(1); \
    } \
} while(0)

// Constants for Carlson functions
__constant__ double ERRTOL_RF = 0.08;
__constant__ double ERRTOL_RC = 0.04; 
__constant__ double ERRTOL_RD = 0.05;
__constant__ double ERRTOL_RJ = 0.05;

__constant__ double TINY = 1.5e-38;
__constant__ double BIG = 3.0e37;

/**
 * Carlson's elliptic integral of the first kind
 * RF(x,y,z) = integral from 0 to infinity of:
 *   (1/2) * (t+x)^(-1/2) * (t+y)^(-1/2) * (t+z)^(-1/2) dt
 * 
 * Where x,y,z >= 0 and at most one of them is zero.
 */
__device__ double carlson_rf(double x, double y, double z) {
    // Constants
    const double third = 1.0/3.0;
    const double c1 = 1.0/24.0;
    const double c2 = 0.1;
    const double c3 = 3.0/44.0;
    const double c4 = 1.0/14.0;
    
    // Validity check (in production, may want to handle gracefully)
    if (fmin(fmin(x,y),z) < 0.0 || 
        fmin(fmin(x+y, x+z), y+z) < TINY ||
        fmax(fmax(x,y),z) > BIG) {
        return NAN; // Invalid arguments
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
        if (fabs(ave) < 1e-15) return NAN;  // Avoid division by zero
        
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        
        if (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL_RF) {
            double e2 = delx * dely - delz * delz;
            double e3 = delx * dely * delz;
            return (1.0 + (c1*e2 - c2 - c3*e3)*e2 + c4*e3) / sqrt(ave);
        }
    }
    
    // Convergence failed after max iterations
    return NAN;
}

/**
 * Carlson's degenerate elliptic integral
 * RC(x,y) = integral from 0 to infinity of:
 *   (1/2) * (t+x)^(-1/2) * (t+y)^(-1) dt
 * 
 * Where x >= 0 and y != 0.
 */
__device__ double carlson_rc(double x, double y) {
    // Constants
    const double third = 1.0/3.0;
    const double c1 = 0.3;
    const double c2 = 1.0/7.0;
    const double c3 = 0.375;
    const double c4 = 9.0/22.0;
    const double sqrtny = 1.3e-19;
    const double tnbg = TINY * BIG;
    const double comp1 = 2.236 / sqrtny;
    const double comp2 = tnbg * tnbg / 25.0;
    
    // Validity check
    if (x < 0.0 || y == 0.0 || 
        (x + fabs(y)) < TINY ||
        (x + fabs(y)) > BIG ||
        (y < -comp1 && x > 0.0 && x < comp2)) {
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
        if (fabs(ave) < 1e-15) return NAN;  // Avoid division by zero
        
        double s = (yt - ave) / ave;
        
        if (fabs(s) <= ERRTOL_RC) {
            return w * (1.0 + s*s*(c1 + s*(c2 + s*(c3 + s*c4)))) / sqrt(ave);
        }
    }
    
    // Convergence failed after max iterations
    return NAN;
}

/**
 * Carlson's elliptic integral of the second kind
 * RD(x,y,z) = integral from 0 to infinity of:
 *   (3/2) * (t+x)^(-1/2) * (t+y)^(-1/2) * (t+z)^(-3/2) dt
 * 
 * Where x,y >= 0, at most one of them is zero, and z > 0.
 */
__device__ double carlson_rd(double x, double y, double z) {
    // Constants  
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/6.0;
    const double c3 = 9.0/22.0;
    const double c4 = 3.0/26.0;
    const double c5 = 0.25 * c3;
    const double c6 = 1.5 * c4;
    const double tiny_rd = 1.0e-25;
    const double big_rd = 4.5e21;
    
    // Validity check
    if (fmin(x,y) < 0.0 || 
        fmin(x+y, z) < tiny_rd ||
        fmax(fmax(x,y),z) > big_rd) {
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
        
        double denominator = sqrtz * (zt + alamb);
        if (fabs(denominator) > 1e-15) {
            sum = sum + fac / denominator;
        }
        fac = 0.25 * fac;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        
        double ave = 0.2 * (xt + yt + 3.0 * zt);
        if (fabs(ave) < 1e-15) return NAN;  // Avoid division by zero
        
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        
        if (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL_RD) {
            double ea = delx * dely;
            double eb = delz * delz;
            double ec = ea - eb;
            double ed = ea - 6.0 * eb;
            double ee = ed + ec + ec;
            
            return 3.0 * sum + fac * (1.0 + ed*(-c1 + c5*ed - c6*delz*ee) +
                   delz*(c2*ee + delz*(-c3*ec + delz*c4*ea))) / (ave * sqrt(ave));
        }
    }
    
    // Convergence failed after max iterations
    return NAN;
}

/**
 * Carlson's elliptic integral of the third kind
 * RJ(x,y,z,p) = integral from 0 to infinity of:
 *   (3/2) * (t+x)^(-1/2) * (t+y)^(-1/2) * (t+z)^(-1/2) * (t+p)^(-1) dt
 * 
 * Where x,y,z >= 0, at most one of them is zero, and p > 0.
 */
__device__ double carlson_rj(double x, double y, double z, double p) {
    // Constants
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/3.0;
    const double c3 = 3.0/22.0;
    const double c4 = 3.0/26.0;
    const double c5 = 0.75 * c3;
    const double c6 = 1.5 * c4;
    const double c7 = 0.5 * c2;
    const double c8 = c3 + c3;
    const double tiny_rj = 2.5e-13;
    const double big_rj = 9.0e11;
    
    // Validity check
    if (fmin(fmin(x,y),z) < 0.0 || 
        fmin(fmin(fmin(x+y, x+z), y+z), p) < tiny_rj ||
        fmax(fmax(fmax(x,y),z),p) > big_rj) {
        return NAN;
    }
    
    double sum = 0.0;
    double fac = 1.0;
    double xt, yt, zt, pt;
    double rcx = 0.0;
    
    if (p > 0.0) {
        xt = x;
        yt = y; 
        zt = z;
        pt = p;
    } else {
        xt = fmin(fmin(x,y),z);
        zt = fmax(fmax(x,y),z);
        yt = x + y + z - xt - zt;
        double a = 1.0 / (yt - p);
        double b = a * (zt - yt) * (yt - xt);
        pt = yt + b;
        double rho = xt * zt / yt;
        double tau = p * pt / yt;
        rcx = carlson_rc(rho, tau);
    }
    
    for (int iter = 0; iter < 100; iter++) {
        double sqrtx = sqrt(xt);
        double sqrty = sqrt(yt);
        double sqrtz = sqrt(zt);
        double sqrtp = sqrt(pt);
        double dnm = sqrtp * (sqrtp + sqrtx) * (sqrtp + sqrty) * (sqrtp + sqrtz);
        
        if (fabs(dnm) > 1e-15) {
            sum = sum + fac / dnm;
        }
        fac = 0.25 * fac;
        double alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        double alpha = (pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
        alpha = alpha * alpha;
        double beta = pt * (pt + alamb);
        beta = beta * beta;
        
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        
        double ave = 0.2 * (xt + yt + zt + pt + pt);
        if (fabs(ave) < 1e-15) return NAN;  // Avoid division by zero
        
        double delx = (ave - xt) / ave;
        double dely = (ave - yt) / ave;
        double delz = (ave - zt) / ave;
        double delp = (ave - pt) / ave;
        
        if (fmax(fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)), fabs(delp)) <= ERRTOL_RJ) {
            double ea = delx * (dely + delz) + dely * delz;
            double eb = delx * dely * delz;
            double ec = delp * delp;
            double ed = ea - 3.0 * ec;
            double ee = eb + 2.0 * delp * (ea - ec);
            
            double result = 3.0 * sum + fac * (1.0 + ed*(-c1 + c5*ed - c6*ee) + 
                           eb*(c7 + delp*(-c8 + delp*c4)) + 
                           delp*ea*(c2 - delp*c3) - c2*delp*ec) / (ave * sqrt(ave));
            
            if (p <= 0.0) {
                double a = 1.0 / (yt - p);
                double b = a * (zt - yt) * (yt - xt);
                result = a * (b * result + 3.0 * (rcx - carlson_rf(xt, yt, zt)));
            }
            return result;
        }
    }
    
    // Convergence failed after max iterations
    return NAN;
}

/**
 * CUDA kernel to batch process elliptic integral test cases
 */
__global__ void elliptic_batch_kernel(
    const char* test_types, 
    const double* x_vals, const double* y_vals, 
    const double* z_vals, const double* p_vals,
    double* results, int n_tests
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_tests) return;
    
    // Extract test type (2 characters)
    char type[3];
    type[0] = test_types[idx * 2];
    type[1] = test_types[idx * 2 + 1];
    type[2] = '\0';
    
    double x = x_vals[idx];
    double y = y_vals[idx];
    double z = z_vals[idx];
    double p = p_vals[idx];
    
    if (type[0] == 'R' && type[1] == 'F') {
        results[idx] = carlson_rf(x, y, z);
    } else if (type[0] == 'R' && type[1] == 'C') {
        results[idx] = carlson_rc(x, y);
    } else if (type[0] == 'R' && type[1] == 'D') {
        results[idx] = carlson_rd(x, y, z);
    } else if (type[0] == 'R' && type[1] == 'J') {
        results[idx] = carlson_rj(x, y, z, p);
    } else {
        results[idx] = NAN; // Unknown test type
    }
}

/**
 * Host function to validate CUDA implementation against FORTRAN reference
 */
extern "C" void validate_carlson_cuda() {
    std::cout << "CUDA Carlson Elliptic Integral Validation" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // This would load test cases and reference data, run CUDA kernels,
    // and compare results. Implementation depends on data loading utilities.
    
    std::cout << "Validation stub - implement data loading and comparison" << std::endl;
}

// Simple test function for development
extern "C" void test_carlson_cuda() {
    // Test RF(1,2,0) = 1.3110287771461 (approximate)
    double x = 1.0, y = 2.0, z = 0.0;
    
    double *d_x, *d_y, *d_z, *d_p, *d_result;
    char *d_types;
    
    CUDA_CHECK(cudaMalloc(&d_x, sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_y, sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_z, sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_p, sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_result, sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_types, 2));
    
    CUDA_CHECK(cudaMemcpy(d_x, &x, sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_y, &y, sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z, &z, sizeof(double), cudaMemcpyHostToDevice));
    
    double p = 0.0;
    char types[2] = {'R', 'F'};
    CUDA_CHECK(cudaMemcpy(d_p, &p, sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_types, types, 2, cudaMemcpyHostToDevice));
    
    elliptic_batch_kernel<<<1, 1>>>(d_types, d_x, d_y, d_z, d_p, d_result, 1);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    double result;
    CUDA_CHECK(cudaMemcpy(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost));
    
    std::cout << "RF(1, 2, 0) = " << result << " (expected ~1.311)" << std::endl;
    
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(d_p);
    cudaFree(d_result);
    cudaFree(d_types);
}