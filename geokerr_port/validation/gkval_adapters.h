/*
 * GEOKERR AND CUDA PORT SOLVER ADAPTERS
 * 
 * Implements the solver interface for both the original GeoKerr
 * FORTRAN implementation and the CUDA port, providing a unified
 * streaming API for validation.
 */

#ifndef GKVAL_ADAPTERS_H
#define GKVAL_ADAPTERS_H

#include "gkval_validator.h"

#ifdef __cplusplus
extern "C" {
#endif

// =============================================================================
// GEOKERR FORTRAN ADAPTER
// =============================================================================

typedef struct {
    // FORTRAN executable path
    char executable_path[512];
    
    // Current case state
    gkval_inputs_t current_inputs;
    bool case_initialized;
    
    // Grid parameters
    uint64_t grid_size;
    double* grid_points;
    char grid_type[64];
    double grid_delta;
    // Per-case grid info
    double lam0;
    double dlam;
    int tpmi;
    int tpri;
    
    // Cached field data
    double* field_data[GKVAL_MAX_FIELDS];
    bool field_computed[GKVAL_MAX_FIELDS];
    
    // Working directory for temporary files
    char work_dir[512];
    
} gkval_geokerr_ctx_t;

// Create GeoKerr FORTRAN solver
gkval_solver_t* gkval_create_geokerr_solver(const char* executable_path,
                                           const char* work_dir,
                                           const gkval_grid_t* grid);

// =============================================================================
// CUDA PORT ADAPTER  
// =============================================================================

typedef struct {
    // CUDA context and stream
    void* cuda_context;
    void* cuda_stream;
    
    // Current case state
    gkval_inputs_t current_inputs;
    bool case_initialized;
    
    // Grid parameters
    uint64_t grid_size;
    double* grid_points;
    
    // Device memory buffers
    double* d_field_data[GKVAL_MAX_FIELDS];
    double* h_field_data[GKVAL_MAX_FIELDS];
    bool field_computed[GKVAL_MAX_FIELDS];
    
    // CUDA solver configuration
    int device_id;
    size_t max_segments_per_kernel;
    
} gkval_cuda_ctx_t;

// Create CUDA port solver
gkval_solver_t* gkval_create_cuda_solver(int device_id,
                                        const gkval_grid_t* grid);

// =============================================================================
// GRID MANAGEMENT
// =============================================================================

// Generate uniform grid in lambda (Mino time parameter)
int gkval_generate_uniform_lambda_grid(const gkval_inputs_t* inputs,
                                      uint64_t S, double delta,
                                      double** grid_points);

// Generate uniform grid in affine parameter
int gkval_generate_uniform_affine_grid(const gkval_inputs_t* inputs,
                                      uint64_t S, double delta,
                                      double** grid_points);

// Generate adaptive grid in lambda based on Tr/Ttheta with policy K*max(periods)
int gkval_generate_adaptive_lambda_grid(const gkval_inputs_t* inputs,
                                       uint64_t S, double K,
                                       double** grid_points,
                                       double* out_dlam);

// Resample field data to a different grid (monotone cubic Hermite)
int gkval_resample_field_data(const double* src_grid, const double* src_data, size_t src_n,
                             const double* dst_grid, double* dst_data, size_t dst_n);

// =============================================================================
// SOLVER INTERFACE IMPLEMENTATIONS
// =============================================================================

// GeoKerr FORTRAN solver functions
int gkval_geokerr_init_case(gkval_solver_t* solver, const gkval_inputs_t* inputs);
int gkval_geokerr_read_field(gkval_solver_t* solver, gkval_field_t field,
                            uint64_t start, uint64_t n, double* output);
void gkval_geokerr_cleanup_case(gkval_solver_t* solver);
void gkval_geokerr_cleanup(gkval_solver_t* solver);

// Accessors for grid and turning points (reference solver)
int gkval_geokerr_get_grid_info(const gkval_solver_t* solver, double* lam0, double* dlam);
int gkval_geokerr_get_tp_indices(const gkval_solver_t* solver, int* tpmi, int* tpri);

// CUDA port solver functions  
int gkval_cuda_init_case(gkval_solver_t* solver, const gkval_inputs_t* inputs);
int gkval_cuda_read_field(gkval_solver_t* solver, gkval_field_t field,
                         uint64_t start, uint64_t n, double* output);
void gkval_cuda_cleanup_case(gkval_solver_t* solver);
void gkval_cuda_cleanup(gkval_solver_t* solver);

// =============================================================================
// UTILITIES
// =============================================================================

// Execute FORTRAN binary with inputs and parse output
int gkval_execute_geokerr_fortran(const char* executable, const char* work_dir,
                                 const gkval_inputs_t* inputs, const double* grid,
                                 uint64_t grid_size, double lam0, double dlam,
                                 int* out_tpmi, int* out_tpri,
                                 double* field_data[GKVAL_MAX_FIELDS]);

// Launch CUDA kernel for geodesic integration
int gkval_launch_cuda_integration(gkval_cuda_ctx_t* ctx, const gkval_inputs_t* inputs,
                                 const double* grid, uint64_t grid_size);

// Convert field data between different coordinate systems
int gkval_transform_coordinates(gkval_field_t from_system, gkval_field_t to_system,
                               const double* input, double* output, size_t n,
                               const gkval_inputs_t* params);

#ifdef __cplusplus
}
#endif

#endif // GKVAL_ADAPTERS_H