/*
 * GEOKERR AND CUDA PORT SOLVER ADAPTERS IMPLEMENTATION
 * 
 * Provides streaming interfaces to both the FORTRAN reference
 * and CUDA port implementations.
 */

#include "gkval_adapters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// Define M_PI if not provided by math.h
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

// Local FNV-1a 64-bit hash for cache keys (file-scope)
static inline uint64_t fnv1a64(const void* data, size_t len) {
    const uint8_t* p = (const uint8_t*)data;
    uint64_t h = 1469598103934665603ULL; // FNV offset basis
    for (size_t i = 0; i < len; ++i) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ULL; // FNV prime
    }
    return h;
}

// ISO_C_BINDING FORTRAN shared-lib entry points
extern void geokerr_eval_lambda(
    double u0, double mu0, double uf, double a, double Lz, double Q2, double t0,
    double lam0, double dlam, int S,
    double* u, double* mu, double* t, double* phi, double* affine,
    int* tpmi, int* tpri, int* info);

extern void geokerr_estimate_periods(
    double u0, double mu0, double a, double Lz, double Q2,
    double* Tr, double* Ttheta, int* info);

extern double geokerr_phi_omega(
    double u0, double mu0, double a, double Lz, double Q2);

// =============================================================================
// GRID GENERATION
// =============================================================================

int gkval_generate_uniform_lambda_grid(const gkval_inputs_t* inputs,
                                      uint64_t S, double delta,
                                      double** grid_points) {
    // inputs is unused for uniform grid; allow NULL
    if (!grid_points || S == 0) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    *grid_points = malloc(S * sizeof(double));
    if (!*grid_points) {
        return GKVAL_ERROR_MEMORY;
    }
    
    // Generate uniform grid in Mino time parameter λ
    // λ spans from 0 to some maximum value based on the geodesic
    double lambda_max = delta * (S - 1);
    (void)lambda_max;
    
    for (uint64_t i = 0; i < S; i++) {
        (*grid_points)[i] = i * delta;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_generate_uniform_affine_grid(const gkval_inputs_t* inputs,
                                      uint64_t S, double delta,
                                      double** grid_points) {
    // inputs is unused for uniform affine grid; allow NULL
    if (!grid_points || S == 0) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    *grid_points = malloc(S * sizeof(double));
    if (!*grid_points) {
        return GKVAL_ERROR_MEMORY;
    }
    
    // Generate uniform grid in affine parameter
    for (uint64_t i = 0; i < S; i++) {
        (*grid_points)[i] = i * delta;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_generate_adaptive_lambda_grid(const gkval_inputs_t* inputs,
                                       uint64_t S, double K,
                                       double** grid_points,
                                       double* out_dlam) {
    if (!inputs || !grid_points || S == 0 || !out_dlam) return GKVAL_ERROR_INVALID_ARG;

    double Tr = -1.0, Tth = -1.0;
    int info = 0;
    geokerr_estimate_periods(inputs->u0, inputs->mu0, inputs->a, inputs->lam, inputs->q2,
                             &Tr, &Tth, &info);
    if (info != 0) {
        // Fallback: simple default
        Tr = 1.0; Tth = 1.0;
    }

    // Degeneracy rules
    const double eps = 1e-14;
    double Tchar;
    if (fabs(inputs->q2) < eps && Tr > 0) {
        Tchar = Tr; // equatorial
    } else if (fabs(inputs->uf - inputs->u0) < eps && Tth > 0) {
        Tchar = Tth; // circular-ish
    } else if ((fabs(inputs->q2) < eps) && (fabs(inputs->uf - inputs->u0) < eps)) {
        // Circular equatorial: use K*(2π/Ωφ)
        double omega = geokerr_phi_omega(inputs->u0, inputs->mu0, inputs->a, inputs->lam, inputs->q2);
        if (fabs(omega) < eps) omega = 1.0; // guard
        Tchar = (2.0 * M_PI) / fabs(omega);
    } else {
        Tchar = fmax(Tr, Tth);
    }
    if (Tchar <= 0.0) Tchar = 1.0;

    double lambda_max = K * Tchar;
    double dlam = lambda_max / (double)(S - 1);
    // Round naturally by storing back
    lambda_max = dlam * (double)(S - 1);

    *grid_points = (double*)malloc(S * sizeof(double));
    if (!*grid_points) return GKVAL_ERROR_MEMORY;
    for (uint64_t i = 0; i < S; i++) {
        (*grid_points)[i] = (double)i * dlam;
    }
    *out_dlam = dlam;
    return GKVAL_SUCCESS;
}

// =============================================================================
// GEOKERR FORTRAN ADAPTER
// =============================================================================

gkval_solver_t* gkval_create_geokerr_solver(const char* executable_path,
                                           const char* work_dir,
                                           const gkval_grid_t* grid) {
    if (!executable_path || !work_dir || !grid) {
        return NULL;
    }
    
    gkval_solver_t* solver = malloc(sizeof(gkval_solver_t));
    if (!solver) return NULL;
    
    gkval_geokerr_ctx_t* ctx = malloc(sizeof(gkval_geokerr_ctx_t));
    if (!ctx) {
        free(solver);
        return NULL;
    }
    
    memset(ctx, 0, sizeof(*ctx));
    strncpy(ctx->executable_path, executable_path, sizeof(ctx->executable_path) - 1);
    strncpy(ctx->work_dir, work_dir, sizeof(ctx->work_dir) - 1);
    ctx->grid_size = grid->S;
    strncpy(ctx->grid_type, grid->type, sizeof(ctx->grid_type) - 1);
    ctx->grid_delta = grid->delta;
    // Allocate grid array now; fill per case in init_case
    ctx->grid_points = (double*)malloc(grid->S * sizeof(double));
    if (!ctx->grid_points) { free(ctx); free(solver); return NULL; }
    
    // Allocate field data arrays
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->field_data[i] = malloc(grid->S * sizeof(double));
        if (!ctx->field_data[i]) {
            // Cleanup on failure
            for (int j = 0; j < i; j++) {
                free(ctx->field_data[j]);
            }
            free(ctx->grid_points);
            free(ctx);
            free(solver);
            return NULL;
        }
        ctx->field_computed[i] = false;
    }
    
    // Set up solver interface
    solver->ctx = ctx;
    solver->init_case = gkval_geokerr_init_case;
    solver->read_field = gkval_geokerr_read_field;
    solver->cleanup_case = gkval_geokerr_cleanup_case;
    solver->cleanup = gkval_geokerr_cleanup;
    
    return solver;
}

int gkval_geokerr_init_case(gkval_solver_t* solver, const gkval_inputs_t* inputs) {
    if (!solver || !solver->ctx || !inputs) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    gkval_geokerr_ctx_t* ctx = (gkval_geokerr_ctx_t*)solver->ctx;
    
    // Store inputs
    ctx->current_inputs = *inputs;
    ctx->case_initialized = true;
    
    // Generate per-case grid
    int gres = GKVAL_SUCCESS;
    if (strcmp(ctx->grid_type, "adaptive_lambda") == 0) {
        double dlam = 0.0;
        double* tmp = NULL;
        gres = gkval_generate_adaptive_lambda_grid(inputs, ctx->grid_size, 6.0, &tmp, &dlam);
        if (gres == GKVAL_SUCCESS) {
            memcpy(ctx->grid_points, tmp, ctx->grid_size * sizeof(double));
            free(tmp);
            ctx->lam0 = (ctx->grid_size > 0) ? ctx->grid_points[0] : 0.0;
            ctx->dlam = dlam;
        }
    } else if (strcmp(ctx->grid_type, "uniform_lambda") == 0) {
        double* tmp = NULL;
        gres = gkval_generate_uniform_lambda_grid(inputs, ctx->grid_size, ctx->grid_delta, &tmp);
        if (gres == GKVAL_SUCCESS) {
            memcpy(ctx->grid_points, tmp, ctx->grid_size * sizeof(double));
            free(tmp);
            ctx->lam0 = (ctx->grid_size > 0) ? ctx->grid_points[0] : 0.0;
            ctx->dlam = (ctx->grid_size > 1) ? (ctx->grid_points[1] - ctx->grid_points[0]) : 0.0;
        }
    } else if (strcmp(ctx->grid_type, "uniform_affine") == 0) {
        double* tmp = NULL;
        gres = gkval_generate_uniform_affine_grid(inputs, ctx->grid_size, ctx->grid_delta, &tmp);
        if (gres == GKVAL_SUCCESS) {
            memcpy(ctx->grid_points, tmp, ctx->grid_size * sizeof(double));
            free(tmp);
            ctx->lam0 = (ctx->grid_size > 0) ? ctx->grid_points[0] : 0.0;
            ctx->dlam = (ctx->grid_size > 1) ? (ctx->grid_points[1] - ctx->grid_points[0]) : 0.0;
        }
    } else {
        return GKVAL_ERROR_INVALID_ARG;
    }
    if (gres != GKVAL_SUCCESS) return gres;

    // Mark all fields as not computed
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->field_computed[i] = false;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_geokerr_read_field(gkval_solver_t* solver, gkval_field_t field,
                            uint64_t start, uint64_t n, double* output) {
    if (!solver || !solver->ctx || !output) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    gkval_geokerr_ctx_t* ctx = (gkval_geokerr_ctx_t*)solver->ctx;
    
    if (!ctx->case_initialized || field >= GKVAL_MAX_FIELDS) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    if (start + n > ctx->grid_size) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    // Compute field data if not already done
    if (!ctx->field_computed[field]) {
        int tpmi = -1, tpri = -1;
        int ret = gkval_execute_geokerr_fortran(ctx->executable_path, ctx->work_dir,
                                               &ctx->current_inputs, ctx->grid_points,
                                               ctx->grid_size, ctx->lam0, ctx->dlam,
                                               &tpmi, &tpri,
                                               ctx->field_data);
        if (ret != GKVAL_SUCCESS) {
            return ret;
        }
        ctx->tpmi = tpmi;
        ctx->tpri = tpri;
        
        // Mark all fields as computed (FORTRAN computes everything at once)
        for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
            ctx->field_computed[i] = true;
        }
    }
    
    // Copy requested segment
    memcpy(output, ctx->field_data[field] + start, n * sizeof(double));
    
    return GKVAL_SUCCESS;
}

void gkval_geokerr_cleanup_case(gkval_solver_t* solver) {
    if (!solver || !solver->ctx) return;
    
    gkval_geokerr_ctx_t* ctx = (gkval_geokerr_ctx_t*)solver->ctx;
    ctx->case_initialized = false;
    
    // Mark fields as not computed for next case
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->field_computed[i] = false;
    }
}

void gkval_geokerr_cleanup(gkval_solver_t* solver) {
    if (!solver || !solver->ctx) return;
    
    gkval_geokerr_ctx_t* ctx = (gkval_geokerr_ctx_t*)solver->ctx;
    
    // Free field data arrays
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        free(ctx->field_data[i]);
    }
    free(ctx->grid_points);
    free(ctx);
    free(solver);
}

int gkval_geokerr_get_grid_info(const gkval_solver_t* solver, double* lam0, double* dlam) {
    if (!solver || !solver->ctx || !lam0 || !dlam) return GKVAL_ERROR_INVALID_ARG;
    const gkval_geokerr_ctx_t* ctx = (const gkval_geokerr_ctx_t*)solver->ctx;
    *lam0 = ctx->lam0;
    *dlam = ctx->dlam;
    return GKVAL_SUCCESS;
}

int gkval_geokerr_get_tp_indices(const gkval_solver_t* solver, int* tpmi, int* tpri) {
    if (!solver || !solver->ctx || !tpmi || !tpri) return GKVAL_ERROR_INVALID_ARG;
    const gkval_geokerr_ctx_t* ctx = (const gkval_geokerr_ctx_t*)solver->ctx;
    *tpmi = ctx->tpmi;
    *tpri = ctx->tpri;
    return GKVAL_SUCCESS;
}

// =============================================================================
// CUDA PORT ADAPTER
// =============================================================================

gkval_solver_t* gkval_create_cuda_solver(int device_id, const gkval_grid_t* grid) {
    if (!grid) return NULL;
    
    gkval_solver_t* solver = malloc(sizeof(gkval_solver_t));
    if (!solver) return NULL;
    
    gkval_cuda_ctx_t* ctx = malloc(sizeof(gkval_cuda_ctx_t));
    if (!ctx) {
        free(solver);
        return NULL;
    }
    
    memset(ctx, 0, sizeof(*ctx));
    ctx->device_id = device_id;
    ctx->grid_size = grid->S;
    ctx->max_segments_per_kernel = 1024; // Default
    
    // Generate grid points (adaptive falls back to uniform spacing for now)
    int gres = GKVAL_SUCCESS;
    if (strcmp(grid->type, "uniform_lambda") == 0) {
        gres = gkval_generate_uniform_lambda_grid(NULL, grid->S, grid->delta, &ctx->grid_points);
    } else if (strcmp(grid->type, "adaptive_lambda") == 0) {
        gres = gkval_generate_uniform_lambda_grid(NULL, grid->S, grid->delta, &ctx->grid_points);
    } else if (strcmp(grid->type, "uniform_affine") == 0) {
        gres = gkval_generate_uniform_affine_grid(NULL, grid->S, grid->delta, &ctx->grid_points);
    } else {
        free(ctx);
        free(solver);
        return NULL;
    }
    if (gres != GKVAL_SUCCESS || ctx->grid_points == NULL) {
        free(ctx);
        free(solver);
        return NULL;
    }
    
    // Allocate host field data arrays
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->h_field_data[i] = malloc(grid->S * sizeof(double));
        if (!ctx->h_field_data[i]) {
            // Cleanup on failure
            for (int j = 0; j < i; j++) {
                free(ctx->h_field_data[j]);
            }
            free(ctx->grid_points);
            free(ctx);
            free(solver);
            return NULL;
        }
        ctx->field_computed[i] = false;
    }
    
    // TODO: Initialize CUDA context and allocate device memory
    
    // Set up solver interface
    solver->ctx = ctx;
    solver->init_case = gkval_cuda_init_case;
    solver->read_field = gkval_cuda_read_field;
    solver->cleanup_case = gkval_cuda_cleanup_case;
    solver->cleanup = gkval_cuda_cleanup;
    
    return solver;
}

int gkval_cuda_init_case(gkval_solver_t* solver, const gkval_inputs_t* inputs) {
    if (!solver || !solver->ctx || !inputs) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    gkval_cuda_ctx_t* ctx = (gkval_cuda_ctx_t*)solver->ctx;
    
    // Store inputs
    ctx->current_inputs = *inputs;
    ctx->case_initialized = true;
    
    // Regenerate per-case grid matching the adaptive policy (K=6) to align with reference
    if (ctx->grid_points && ctx->grid_size > 0) {
        double* tmp = NULL;
        double dlam = 0.0;
        int gres = gkval_generate_adaptive_lambda_grid(inputs, ctx->grid_size, 6.0, &tmp, &dlam);
        if (gres == GKVAL_SUCCESS && tmp) {
            memcpy(ctx->grid_points, tmp, ctx->grid_size * sizeof(double));
            free(tmp);
        } else {
            // Fallback to uniform if adaptive fails
            if (ctx->grid_size > 0) {
                for (uint64_t i = 0; i < ctx->grid_size; i++) ctx->grid_points[i] = (double)i * 1e-3;
            }
        }
    }

    // Mark all fields as not computed
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->field_computed[i] = false;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_cuda_read_field(gkval_solver_t* solver, gkval_field_t field,
                         uint64_t start, uint64_t n, double* output) {
    if (!solver || !solver->ctx || !output) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    gkval_cuda_ctx_t* ctx = (gkval_cuda_ctx_t*)solver->ctx;
    
    if (!ctx->case_initialized || field >= GKVAL_MAX_FIELDS) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    if (start + n > ctx->grid_size) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    // Compute field data if not already done
    if (!ctx->field_computed[field]) {
        int ret = gkval_launch_cuda_integration(ctx, &ctx->current_inputs,
                                               ctx->grid_points, ctx->grid_size);
        if (ret != GKVAL_SUCCESS) {
            return ret;
        }
        
        // Mark all fields as computed (CUDA computes everything at once)
        for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
            ctx->field_computed[i] = true;
        }
    }
    
    // Copy requested segment
    memcpy(output, ctx->h_field_data[field] + start, n * sizeof(double));
    
    return GKVAL_SUCCESS;
}

void gkval_cuda_cleanup_case(gkval_solver_t* solver) {
    if (!solver || !solver->ctx) return;
    
    gkval_cuda_ctx_t* ctx = (gkval_cuda_ctx_t*)solver->ctx;
    ctx->case_initialized = false;
    
    // Mark fields as not computed for next case
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        ctx->field_computed[i] = false;
    }
}

void gkval_cuda_cleanup(gkval_solver_t* solver) {
    if (!solver || !solver->ctx) return;
    
    gkval_cuda_ctx_t* ctx = (gkval_cuda_ctx_t*)solver->ctx;
    
    // Free host field data arrays
    for (int i = 0; i < GKVAL_MAX_FIELDS; i++) {
        free(ctx->h_field_data[i]);
    }
    
    // TODO: Free CUDA device memory and context
    
    free(ctx->grid_points);
    free(ctx);
    free(solver);
}

// =============================================================================
// EXECUTION FUNCTIONS (PLACEHOLDERS)
// =============================================================================

int gkval_execute_geokerr_fortran(const char* executable, const char* work_dir,
                                 const gkval_inputs_t* inputs, const double* grid,
                                 uint64_t grid_size, double lam0, double dlam,
                                 int* out_tpmi, int* out_tpri,
                                 double* field_data[GKVAL_MAX_FIELDS]) {
    if (!inputs || !grid || !field_data || grid_size == 0) {
        return GKVAL_ERROR_INVALID_ARG;
    }

    // Simple on-disk cache to avoid recomputation across runs for the
    // first up to 10k samples (covers default S=4096).
    // Cache key includes inputs, lam0, dlam, and S.
    enum { GKVAL_GEOKERR_CACHE_K = 10000 };

    // Build a compact key struct
    struct key_t {
        gkval_inputs_t in;
        double lam0;
        double dlam;
        uint64_t S;
    } key;
    key.in = *inputs;
    key.lam0 = lam0;
    key.dlam = dlam;
    key.S = grid_size;
    uint64_t h = fnv1a64(&key, sizeof(key));

    // Resolve cache path: ./cache/geokerr/<hexhash>.bin
    char cache_dir[256] = {0};
    snprintf(cache_dir, sizeof(cache_dir), "cache/geokerr");
    // Best-effort: create directories
    (void)mkdir("cache", 0755);
    (void)mkdir(cache_dir, 0755);

    char cache_path[512];
    snprintf(cache_path, sizeof(cache_path), "%s/%016llx.bin", cache_dir, (unsigned long long)h);

    // Cache header
    typedef struct {
        uint32_t magic;   // 'GKCH'
        uint32_t version; // 1
        uint64_t S;       // number of samples stored
        double lam0;
        double dlam;
        gkval_inputs_t in;
        int32_t tpmi;
        int32_t tpri;
    } cache_header_t;

    // Try load from cache if S fits within cache limit
    if (grid_size <= GKVAL_GEOKERR_CACHE_K) {
        FILE* fc = fopen(cache_path, "rb");
        if (fc) {
            cache_header_t hdr;
            size_t ok = fread(&hdr, sizeof(hdr), 1, fc);
            if (ok == 1 && hdr.magic == 0x474B4348U /*GKCH*/ && hdr.version == 1 && hdr.S == grid_size) {
                // Basic sanity on parameters
                if (hdr.lam0 == lam0 && hdr.dlam == dlam &&
                    memcmp(&hdr.in, inputs, sizeof(gkval_inputs_t)) == 0) {
                    // Read fields in order u, mu, t, phi, affine
                    int good = 1;
                    for (int f = 0; f < GKVAL_MAX_FIELDS; ++f) {
                        if (fread(field_data[f], sizeof(double), grid_size, fc) != grid_size) { good = 0; break; }
                    }
                    if (good) {
                        if (out_tpmi) *out_tpmi = hdr.tpmi;
                        if (out_tpri) *out_tpri = hdr.tpri;
                        fclose(fc);
                        return GKVAL_SUCCESS;
                    }
                }
            }
            fclose(fc);
            // Fall-through to compute and overwrite cache
        }
    }

    // Use shared library entrypoint for per-λ sampling
    // Note: gkval_inputs_t.lam is interpreted as Lz (angular momentum)
    (void)grid; // grid not needed beyond lam0/dlam

    // Temporary buffers
    double* bu = (double*)malloc(grid_size * sizeof(double));
    double* bmu = (double*)malloc(grid_size * sizeof(double));
    double* bt = (double*)malloc(grid_size * sizeof(double));
    double* bphi = (double*)malloc(grid_size * sizeof(double));
    double* baff = (double*)malloc(grid_size * sizeof(double));
    if (!bu || !bmu || !bt || !bphi || !baff) {
        free(bu); free(bmu); free(bt); free(bphi); free(baff);
        return GKVAL_ERROR_MEMORY;
    }

    int tpmi = -1, tpri = -1, info = -1;
    geokerr_eval_lambda(
        inputs->u0, inputs->mu0, inputs->uf, inputs->a, inputs->lam, inputs->q2, inputs->t0,
        lam0, dlam, (int)grid_size,
        bu, bmu, bt, bphi, baff,
        &tpmi, &tpri, &info);

    if (out_tpmi) *out_tpmi = tpmi;
    if (out_tpri) *out_tpri = tpri;

    if (info != 0) {
        // Fill with NaN on failure
        for (uint64_t i = 0; i < grid_size; i++) {
            field_data[GKVAL_FIELD_U][i] = NAN;
            field_data[GKVAL_FIELD_MU][i] = NAN;
            field_data[GKVAL_FIELD_T][i] = NAN;
            field_data[GKVAL_FIELD_PHI][i] = NAN;
            field_data[GKVAL_FIELD_AFFINE][i] = NAN;
        }
        free(bu); free(bmu); free(bt); free(bphi); free(baff);
        return GKVAL_ERROR_IO;
    }

    // Copy into output fields (phi already unwrapped in Fortran routine)
    memcpy(field_data[GKVAL_FIELD_U], bu,   grid_size * sizeof(double));
    memcpy(field_data[GKVAL_FIELD_MU], bmu, grid_size * sizeof(double));
    memcpy(field_data[GKVAL_FIELD_T], bt,   grid_size * sizeof(double));
    memcpy(field_data[GKVAL_FIELD_PHI], bphi, grid_size * sizeof(double));
    memcpy(field_data[GKVAL_FIELD_AFFINE], baff, grid_size * sizeof(double));

    // Save to cache if S within limit
    if (grid_size <= GKVAL_GEOKERR_CACHE_K) {
        FILE* fcw = fopen(cache_path, "wb");
        if (fcw) {
            cache_header_t hdr = {0};
            hdr.magic = 0x474B4348U; // 'GKCH'
            hdr.version = 1;
            hdr.S = grid_size;
            hdr.lam0 = lam0;
            hdr.dlam = dlam;
            hdr.in = *inputs;
            hdr.tpmi = tpmi;
            hdr.tpri = tpri;
            if (fwrite(&hdr, sizeof(hdr), 1, fcw) == 1) {
                // Write fields in fixed order
                (void)fwrite(bu,   sizeof(double), grid_size, fcw);
                (void)fwrite(bmu,  sizeof(double), grid_size, fcw);
                (void)fwrite(bt,   sizeof(double), grid_size, fcw);
                (void)fwrite(bphi, sizeof(double), grid_size, fcw);
                (void)fwrite(baff, sizeof(double), grid_size, fcw);
                fflush(fcw);
            }
            fclose(fcw);
        }
    }

    free(bu); free(bmu); free(bt); free(bphi); free(baff);
    return GKVAL_SUCCESS;
}

// External CUDA functions
extern int geokerr_complete_stream_cuda(
    double u0, double uf, double mu0, double a, double L, double Q2,
    double lam0, double dlam, int S,
    double* h_u, double* h_mu, double* h_t, double* h_phi, double* h_affine);

int gkval_launch_cuda_integration(gkval_cuda_ctx_t* ctx, const gkval_inputs_t* inputs,
                                 const double* grid, uint64_t grid_size) {
    if (!ctx || !inputs || !grid) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    // Use CUDA geodesic stream integration with semi-analytic method for scientific accuracy
    double lam0 = (grid_size > 0) ? grid[0] : 0.0;
    double dlam = (grid_size > 1) ? (grid[1] - grid[0]) : 0.0;

    int rc = geokerr_complete_stream_cuda(
        inputs->u0, inputs->uf, inputs->mu0, inputs->a, inputs->lam, inputs->q2,
        lam0, dlam, (int)grid_size,
        ctx->h_field_data[GKVAL_FIELD_U],
        ctx->h_field_data[GKVAL_FIELD_MU],
        ctx->h_field_data[GKVAL_FIELD_T],
        ctx->h_field_data[GKVAL_FIELD_PHI],
        ctx->h_field_data[GKVAL_FIELD_AFFINE]);
    if (rc != 0) {
        // Fill with NaN on failure
        for (uint64_t i = 0; i < grid_size; i++) {
            ctx->h_field_data[GKVAL_FIELD_U][i] = NAN;
            ctx->h_field_data[GKVAL_FIELD_MU][i] = NAN;
            ctx->h_field_data[GKVAL_FIELD_T][i] = NAN;
            ctx->h_field_data[GKVAL_FIELD_PHI][i] = NAN;
            ctx->h_field_data[GKVAL_FIELD_AFFINE][i] = NAN;
        }
        return GKVAL_ERROR_IO;
    }
    return GKVAL_SUCCESS;
}

// =============================================================================
// UTILITY FUNCTIONS (PLACEHOLDERS)
// =============================================================================

int gkval_resample_field_data(const double* src_grid, const double* src_data, size_t src_n,
                             const double* dst_grid, double* dst_data, size_t dst_n) {
    // TODO: Implement monotone cubic Hermite interpolation
    return GKVAL_ERROR_INVALID_ARG;
}

int gkval_transform_coordinates(gkval_field_t from_system, gkval_field_t to_system,
                               const double* input, double* output, size_t n,
                               const gkval_inputs_t* params) {
    // TODO: Implement coordinate transformations
    return GKVAL_ERROR_INVALID_ARG;
}