/*
 * GEOKERR VALIDATION CORE LIBRARY
 * 
 * Cryptographic validation framework with ULP-based quantization,
 * Merkle tree proofs, and crash-safe resume capability.
 * 
 * Based on the comprehensive validation specification for
 * bit-stable audit trails with minimal storage.
 */

#ifndef GKVAL_CORE_H
#define GKVAL_CORE_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// =============================================================================
// CONSTANTS AND CONFIGURATION
// =============================================================================

#define GKVAL_VERSION "1.0"
#define GKVAL_HASH_SIZE 32
#define GKVAL_MAX_FIELDS 8
#define GKVAL_MAX_FIELD_NAME 32
#define GKVAL_DEFAULT_SEG_SIZE 1024
#define GKVAL_DEFAULT_QUANT_LSB 8
#define GKVAL_DEFAULT_WINDOW_SIZE 32

// Magic numbers for fail window files
#define GKVAL_FAIL_WINDOW_MAGIC 0x474B5057  // "GKPW"
#define GKVAL_FAIL_WINDOW_VERSION 1

// Field indices (must match manifest.json field order)
typedef enum {
    GKVAL_FIELD_U = 0,
    GKVAL_FIELD_MU = 1, 
    GKVAL_FIELD_T = 2,
    GKVAL_FIELD_PHI = 3,
    GKVAL_FIELD_AFFINE = 4
} gkval_field_t;

// =============================================================================
// DATA STRUCTURES
// =============================================================================

// Hash type (BLAKE3 32-byte output)
typedef struct {
    uint8_t bytes[GKVAL_HASH_SIZE];
} gkval_hash_t;

// Segment leaf for Merkle tree
typedef struct {
    uint64_t start;
    uint64_t end;
    gkval_hash_t root;
} gkval_segment_leaf_t;

// Tolerances for acceptance predicate
typedef struct {
    double atol;           // Absolute tolerance
    double rtol;           // Relative tolerance  
    uint64_t max_ulps;     // Maximum ULP distance
} gkval_tolerances_t;

// Grid specification
typedef struct {
    char type[64];         // "uniform_lambda", "resampled", etc.
    uint64_t S;           // Total grid points
    double delta;         // Grid spacing (if uniform)
} gkval_grid_t;

// Case invariants tracking
typedef struct {
    double max_lambda_drift;
    double max_q2_drift;
    uint64_t tpmi;
    uint64_t tpri;
    uint64_t nan_count;
} gkval_invariants_t;

// Validation manifest
typedef struct {
    char run_id[64];
    char geokerr_commit[64];
    char port_commit[64];
    gkval_grid_t grid;
    gkval_tolerances_t tolerances;
    uint8_t quant_lsb;
    char fields[GKVAL_MAX_FIELDS][GKVAL_MAX_FIELD_NAME];
    size_t num_fields;
} gkval_manifest_t;

// Cursor for crash-safe resume
typedef struct {
    uint64_t case_id;
    gkval_field_t field;
    uint64_t start;
} gkval_cursor_t;

// Case inputs
typedef struct {
    double u0, mu0, uf, a, lam, q2, t0;
} gkval_inputs_t;

// =============================================================================
// CORE FUNCTIONS
// =============================================================================

// Initialize validation library
int gkval_init(void);

// Hash functions with domain separation
gkval_hash_t gkval_hash_leaf(const char* run_id, uint64_t case_id, 
                            gkval_field_t field, uint64_t start, uint64_t end,
                            const uint8_t* quantized_data, size_t data_len);

gkval_hash_t gkval_hash_node(const gkval_hash_t* left, const gkval_hash_t* right);

// ULP-based quantization 
int gkval_quantize_double(double x, uint8_t quant_lsb, uint64_t* quantized);
int gkval_quantize_block(const double* data, size_t n, uint8_t quant_lsb, 
                        uint8_t* output, size_t* output_len);

// ULP distance calculation
uint64_t gkval_ulp_distance(double a, double b);

// Acceptance predicate  
bool gkval_accept_sample(double a, double b, const gkval_tolerances_t* tol);
int gkval_accept_block(const double* a, const double* b, size_t n,
                      const gkval_tolerances_t* tol, size_t* first_fail);

// Angle unwrapping for phi field
double gkval_unwrap_angle(double angle_prev, double angle_curr);
int gkval_unwrap_angles(double* angles, size_t n);

// Merkle tree operations
gkval_hash_t gkval_merkle_from_leaves(const gkval_segment_leaf_t* leaves, size_t n);
gkval_hash_t gkval_merkle_from_hashes(const gkval_hash_t* hashes, size_t n);

// =============================================================================
// FILE I/O AND PERSISTENCE
// =============================================================================

// Manifest operations
int gkval_save_manifest(const char* run_dir, const gkval_manifest_t* manifest);
int gkval_load_manifest(const char* run_dir, gkval_manifest_t* manifest);

// Cursor operations (crash-safe resume)
int gkval_save_cursor(const char* run_dir, const gkval_cursor_t* cursor);
int gkval_load_cursor(const char* run_dir, gkval_cursor_t* cursor);

// Ledger operations (NDJSON append)
int gkval_append_ledger_line(const char* run_dir, const char* json_line);

// Fail window operations
int gkval_save_fail_window(const char* run_dir, uint64_t case_id, 
                          gkval_field_t field, uint64_t start, uint64_t end,
                          const double* geokerr_data, const double* port_data);

int gkval_load_fail_window(const char* run_dir, uint64_t case_id,
                          gkval_field_t field, uint64_t start, uint64_t end,
                          double* geokerr_data, double* port_data);

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

// String conversion utilities
const char* gkval_field_to_string(gkval_field_t field);
gkval_field_t gkval_string_to_field(const char* str);

// Hash utilities
void gkval_hash_to_hex(const gkval_hash_t* hash, char* hex_str);
int gkval_hex_to_hash(const char* hex_str, gkval_hash_t* hash);

// Directory utilities
int gkval_create_run_dir(const char* run_dir);
int gkval_create_fail_window_dir(const char* run_dir);

// Error handling
const char* gkval_strerror(int error_code);

// =============================================================================
// FP ENVIRONMENT (MXCSR) VALIDATION
// =============================================================================

// Ensure x86 SSE control state is set to round-to-nearest and FTZ/DAZ are disabled.
// If fix=true, attempts to set the required flags. Returns GKVAL_SUCCESS if compliant
// or successfully fixed, otherwise GKVAL_ERROR_IO (or GKVAL_SUCCESS on non-x86).
int gkval_validate_fp_env(bool fix);

// =============================================================================
// ERROR CODES
// =============================================================================

#define GKVAL_SUCCESS           0
#define GKVAL_ERROR_INVALID_ARG -1
#define GKVAL_ERROR_IO          -2
#define GKVAL_ERROR_HASH        -3
#define GKVAL_ERROR_QUANTIZE    -4
#define GKVAL_ERROR_MERKLE      -5
#define GKVAL_ERROR_JSON        -6
#define GKVAL_ERROR_MEMORY      -7

#ifdef __cplusplus
}
#endif

#endif // GKVAL_CORE_H