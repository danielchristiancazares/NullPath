/*
 * GEOKERR STREAMING VALIDATOR
 * 
 * Implements the main validation loop with crash-safe resume,
 * progress tracking, and fail window generation.
 */

#ifndef GKVAL_VALIDATOR_H
#define GKVAL_VALIDATOR_H

#include "gkval_core.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// =============================================================================
// STREAMING VALIDATOR DATA STRUCTURES
// =============================================================================

// Solver interface for reading geodesic data
typedef struct gkval_solver gkval_solver_t;

struct gkval_solver {
    void* ctx;  // Solver-specific context
    
    // Initialize solver with case inputs
    int (*init_case)(gkval_solver_t* solver, const gkval_inputs_t* inputs);
    
    // Read field data for a segment [start, start+n)
    int (*read_field)(gkval_solver_t* solver, gkval_field_t field, 
                     uint64_t start, uint64_t n, double* output);
    
    // Cleanup case resources
    void (*cleanup_case)(gkval_solver_t* solver);
    
    // Cleanup solver
    void (*cleanup)(gkval_solver_t* solver);
};

// Validation session state
typedef struct {
    char run_dir[512];
    gkval_manifest_t manifest;
    gkval_cursor_t cursor;
    
    // Solvers
    gkval_solver_t* geokerr_solver;
    gkval_solver_t* port_solver;
    
    // Runtime state
    FILE* ledger_file;
    FILE* progress_file;
    uint64_t processed_cases;
    uint64_t processed_segments;
    uint64_t processed_samples;
    uint64_t failed_cases;
    
    // Performance tracking
    double start_time;
    double last_checkpoint_time;
    
} gkval_session_t;

// Case validation result
typedef struct {
    uint64_t case_id;
    gkval_inputs_t inputs;
    
    // Per-field segment leaves
    gkval_segment_leaf_t* seg_leaves[GKVAL_MAX_FIELDS];
    size_t num_seg_leaves[GKVAL_MAX_FIELDS];
    
    // Field and case roots
    gkval_hash_t field_roots[GKVAL_MAX_FIELDS];
    gkval_hash_t case_root;
    
    // Validation ranges (for partial validation)
    struct {
        uint64_t start, end;
    } validated_ranges[GKVAL_MAX_FIELDS][16];  // Max 16 ranges per field
    size_t num_validated_ranges[GKVAL_MAX_FIELDS];
    
    // Physics invariants
    gkval_invariants_t invariants;
    // Grid bit patterns (reference grid)
    uint64_t lam0_bits;
    uint64_t dlam_bits;
    
    // Status
    bool success;
    gkval_field_t first_fail_field;
    uint64_t first_fail_index;
    
} gkval_case_result_t;

// =============================================================================
// SESSION MANAGEMENT
// =============================================================================

// Initialize validation session
int gkval_init_session(gkval_session_t* session, const char* run_dir,
                      gkval_solver_t* geokerr_solver, gkval_solver_t* port_solver);

// Load existing session for resume
int gkval_load_session(gkval_session_t* session, const char* run_dir,
                      gkval_solver_t* geokerr_solver, gkval_solver_t* port_solver);

// Save session state and cleanup
int gkval_cleanup_session(gkval_session_t* session);

// =============================================================================
// VALIDATION LOOP
// =============================================================================

// Validate a single case
int gkval_validate_case(gkval_session_t* session, uint64_t case_id,
                       const gkval_inputs_t* inputs, gkval_case_result_t* result);

// Validate multiple cases from a list
int gkval_validate_cases(gkval_session_t* session, const gkval_inputs_t* cases,
                        uint64_t num_cases, uint64_t start_case_id);

// Resume validation from cursor
int gkval_resume_validation(gkval_session_t* session);

// =============================================================================
// PROGRESS AND CHECKPOINTING
// =============================================================================

// Update and save cursor (crash-safe)
int gkval_update_cursor(gkval_session_t* session, uint64_t case_id,
                       gkval_field_t field, uint64_t start);

// Checkpoint progress to disk
int gkval_checkpoint_progress(gkval_session_t* session);

// Print validation statistics
void gkval_print_progress(const gkval_session_t* session);

// =============================================================================
// FAIL WINDOW MANAGEMENT
// =============================================================================

// Create fail window around first failure
int gkval_create_fail_window(gkval_session_t* session, uint64_t case_id,
                            gkval_field_t field, uint64_t fail_index,
                            const double* geokerr_data, const double* port_data,
                            uint64_t data_len);

// =============================================================================
// INVARIANTS TRACKING
// =============================================================================

// Initialize invariants for a case
void gkval_init_invariants(gkval_invariants_t* inv);

// Update invariants with field data
int gkval_update_invariants(gkval_invariants_t* inv, gkval_field_t field,
                           const double* geokerr_data, const double* port_data,
                           uint64_t n);

// Check invariants against tolerances
bool gkval_check_invariants(const gkval_invariants_t* inv,
                           const gkval_tolerances_t* tol);

// =============================================================================
// JSON SERIALIZATION
// =============================================================================

// Serialize case result to JSON line
int gkval_case_result_to_json(const gkval_case_result_t* result,
                             const char* run_id, uint64_t grid_S,
                             char* json_line, size_t max_len);

// Serialize manifest to JSON
int gkval_manifest_to_json(const gkval_manifest_t* manifest,
                          char* json_str, size_t max_len);

// Parse manifest from JSON
int gkval_json_to_manifest(const char* json_str, gkval_manifest_t* manifest);

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

// Get current timestamp
double gkval_get_timestamp(void);

// Create run directory structure
int gkval_create_run_structure(const char* run_dir);

// Atomic file operations
int gkval_write_file_atomic(const char* filepath, const char* content);

#ifdef __cplusplus
}
#endif

#endif // GKVAL_VALIDATOR_H