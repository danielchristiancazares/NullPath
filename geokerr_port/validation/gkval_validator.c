/*
 * GEOKERR STREAMING VALIDATOR IMPLEMENTATION
 * 
 * Main validation loop with crash-safe resume, progress tracking,
 * and comprehensive error handling.
 */

#define _GNU_SOURCE
#include "gkval_validator.h"
#include "gkval_adapters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

double gkval_get_timestamp(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

int gkval_create_run_structure(const char* run_dir) {
    if (gkval_create_run_dir(run_dir) != GKVAL_SUCCESS) {
        return GKVAL_ERROR_IO;
    }
    
    if (gkval_create_fail_window_dir(run_dir) != GKVAL_SUCCESS) {
        return GKVAL_ERROR_IO;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_write_file_atomic(const char* filepath, const char* content) {
    char temp_path[1024];
    snprintf(temp_path, sizeof(temp_path), "%s.tmp", filepath);
    
    FILE* f = fopen(temp_path, "w");
    if (!f) return GKVAL_ERROR_IO;
    
    if (fputs(content, f) == EOF) {
        fclose(f);
        unlink(temp_path);
        return GKVAL_ERROR_IO;
    }
    
    if (fflush(f) != 0 || fsync(fileno(f)) != 0) {
        fclose(f);
        unlink(temp_path);
        return GKVAL_ERROR_IO;
    }
    
    fclose(f);
    
    if (rename(temp_path, filepath) != 0) {
        unlink(temp_path);
        return GKVAL_ERROR_IO;
    }
    
    return GKVAL_SUCCESS;
}

// =============================================================================
// SESSION MANAGEMENT
// =============================================================================

int gkval_init_session(gkval_session_t* session, const char* run_dir,
                      gkval_solver_t* geokerr_solver, gkval_solver_t* port_solver) {
    if (!session || !run_dir || !geokerr_solver || !port_solver) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    memset(session, 0, sizeof(*session));
    strncpy(session->run_dir, run_dir, sizeof(session->run_dir) - 1);
    
    session->geokerr_solver = geokerr_solver;
    session->port_solver = port_solver;
    session->start_time = gkval_get_timestamp();
    session->last_checkpoint_time = session->start_time;
    
    // Create run directory structure
    int ret = gkval_create_run_structure(run_dir);
    if (ret != GKVAL_SUCCESS) {
        return ret;
    }
    
    // Open ledger file for append
    char ledger_path[1024];
    snprintf(ledger_path, sizeof(ledger_path), "%s/ledger.ndjson", run_dir);
    session->ledger_file = fopen(ledger_path, "a");
    if (!session->ledger_file) {
        return GKVAL_ERROR_IO;
    }
    
    // Initialize cursor
    session->cursor.case_id = 0;
    session->cursor.field = GKVAL_FIELD_U;
    session->cursor.start = 0;
    
    return GKVAL_SUCCESS;
}

int gkval_load_session(gkval_session_t* session, const char* run_dir,
                      gkval_solver_t* geokerr_solver, gkval_solver_t* port_solver) {
    int ret = gkval_init_session(session, run_dir, geokerr_solver, port_solver);
    if (ret != GKVAL_SUCCESS) {
        return ret;
    }
    
    // Try to load existing cursor
    gkval_cursor_t cursor;
    if (gkval_load_cursor(run_dir, &cursor) == GKVAL_SUCCESS) {
        session->cursor = cursor;
    }
    
    // Try to load existing manifest
    gkval_manifest_t manifest;
    if (gkval_load_manifest(run_dir, &manifest) == GKVAL_SUCCESS) {
        session->manifest = manifest;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_cleanup_session(gkval_session_t* session) {
    if (!session) return GKVAL_ERROR_INVALID_ARG;
    
    // Save final cursor
    gkval_save_cursor(session->run_dir, &session->cursor);
    
    // Close files
    if (session->ledger_file) {
        fflush(session->ledger_file);
        fsync(fileno(session->ledger_file));
        fclose(session->ledger_file);
        session->ledger_file = NULL;
    }
    
    if (session->progress_file) {
        fclose(session->progress_file);
        session->progress_file = NULL;
    }
    
    return GKVAL_SUCCESS;
}

// =============================================================================
// INVARIANTS TRACKING
// =============================================================================

void gkval_init_invariants(gkval_invariants_t* inv) {
    if (!inv) return;
    
    memset(inv, 0, sizeof(*inv));
    inv->max_lambda_drift = 0.0;
    inv->max_q2_drift = 0.0;
}

int gkval_update_invariants(gkval_invariants_t* inv, gkval_field_t field,
                           const double* geokerr_data, const double* port_data,
                           uint64_t n) {
    if (!inv || !geokerr_data || !port_data) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    // Check for NaN/Inf
    for (uint64_t i = 0; i < n; i++) {
        if (isnan(geokerr_data[i]) || isinf(geokerr_data[i]) ||
            isnan(port_data[i]) || isinf(port_data[i])) {
            inv->nan_count++;
        }
    }
    
    // Field-specific invariant checking
    switch (field) {
        case GKVAL_FIELD_AFFINE:
            // Check lambda drift (if this is the affine parameter)
            for (uint64_t i = 0; i < n; i++) {
                double drift = fabs(geokerr_data[i] - port_data[i]);
                if (drift > inv->max_lambda_drift) {
                    inv->max_lambda_drift = drift;
                }
            }
            break;
            
        case GKVAL_FIELD_U:
        case GKVAL_FIELD_MU:
            // Could check coordinate bounds, etc.
            break;
            
        default:
            break;
    }
    
    return GKVAL_SUCCESS;
}

bool gkval_check_invariants(const gkval_invariants_t* inv,
                           const gkval_tolerances_t* tol) {
    if (!inv || !tol) return false;
    
    // Fail if any NaN/Inf detected
    if (inv->nan_count > 0) return false;
    
    // Check lambda drift within tolerance
    if (inv->max_lambda_drift > tol->atol + tol->rtol) {
        return false;
    }
    
    return true;
}

// =============================================================================
// VALIDATION CORE FUNCTIONS
// =============================================================================

int gkval_validate_case(gkval_session_t* session, uint64_t case_id,
                       const gkval_inputs_t* inputs, gkval_case_result_t* result) {
    if (!session || !inputs || !result) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    memset(result, 0, sizeof(*result));
    result->case_id = case_id;
    result->inputs = *inputs;
    result->success = true;
    
    gkval_init_invariants(&result->invariants);
    
    uint64_t S = session->manifest.grid.S;
    uint64_t seg_size = GKVAL_DEFAULT_SEG_SIZE;
    
    // Initialize both solvers for this case
    if (session->geokerr_solver->init_case(session->geokerr_solver, inputs) != 0 ||
        session->port_solver->init_case(session->port_solver, inputs) != 0) {
        return GKVAL_ERROR_IO;
    }
    
    // Process each field
    for (size_t field_idx = 0; field_idx < session->manifest.num_fields; field_idx++) {
        gkval_field_t field = (gkval_field_t)field_idx;
        
        // Skip if already processed (for resume)
        if (case_id == session->cursor.case_id && field < session->cursor.field) {
            continue;
        }
        
        size_t num_segments = 0;
        gkval_segment_leaf_t* segments = malloc((S / seg_size + 1) * sizeof(gkval_segment_leaf_t));
        if (!segments) {
            result->success = false;
            return GKVAL_ERROR_MEMORY;
        }
        
        // Process segments for this field
        for (uint64_t start = 0; start < S; start += seg_size) {
            // Skip if already processed (for resume within field)
            if (case_id == session->cursor.case_id && 
                field == session->cursor.field && 
                start < session->cursor.start) {
                continue;
            }
            
            uint64_t n = (start + seg_size <= S) ? seg_size : (S - start);
            
            // Allocate segment data
            double* geokerr_data = malloc(n * sizeof(double));
            double* port_data = malloc(n * sizeof(double));
            
            if (!geokerr_data || !port_data) {
                free(geokerr_data);
                free(port_data);
                free(segments);
                result->success = false;
                return GKVAL_ERROR_MEMORY;
            }
            
            // Read data from both solvers
            if (session->geokerr_solver->read_field(session->geokerr_solver, field, start, n, geokerr_data) != 0 ||
                session->port_solver->read_field(session->port_solver, field, start, n, port_data) != 0) {
                free(geokerr_data);
                free(port_data);
                free(segments);
                result->success = false;
                return GKVAL_ERROR_IO;
            }
            
            // Special handling for phi field (unwrap angles)
            if (field == GKVAL_FIELD_PHI) {
                gkval_unwrap_angles(geokerr_data, n);
                gkval_unwrap_angles(port_data, n);
            }
            
            // Check acceptance
            size_t first_fail;
            if (gkval_accept_block(geokerr_data, port_data, n, &session->manifest.tolerances, &first_fail) != GKVAL_SUCCESS) {
                // Create fail window (segment-local buffers; convert to absolute indices in helper).
                gkval_create_fail_window(session, case_id, field, start, first_fail,
                                        geokerr_data, port_data, n);
                
                result->success = false;
                result->first_fail_field = field;
                result->first_fail_index = start + first_fail;
                
                free(geokerr_data);
                free(port_data);
                free(segments);
                return GKVAL_SUCCESS;  // Not an error, just a validation failure
            }
            
            // Update invariants
            gkval_update_invariants(&result->invariants, field, geokerr_data, port_data, n);
            
            // Quantize and hash GeoKerr data (the reference)
            size_t quantized_len = n * sizeof(uint64_t);
            uint8_t* quantized = malloc(quantized_len);
            if (!quantized) {
                free(geokerr_data);
                free(port_data);
                free(segments);
                result->success = false;
                return GKVAL_ERROR_MEMORY;
            }
            
            if (gkval_quantize_block(geokerr_data, n, session->manifest.quant_lsb,
                                    quantized, &quantized_len) != GKVAL_SUCCESS) {
                free(geokerr_data);
                free(port_data);
                free(quantized);
                free(segments);
                result->success = false;
                return GKVAL_ERROR_QUANTIZE;
            }
            
            // Create segment leaf hash
            gkval_hash_t leaf_hash = gkval_hash_leaf(session->manifest.run_id, case_id,
                                                    field, start, start + n,
                                                    quantized, quantized_len);
            
            // Store segment
            segments[num_segments].start = start;
            segments[num_segments].end = start + n;
            segments[num_segments].root = leaf_hash;
            num_segments++;
            
            // Update progress
            session->processed_segments++;
            session->processed_samples += n;
            
            // Update cursor
            gkval_update_cursor(session, case_id, field, start + seg_size);
            
            // Cleanup
            free(geokerr_data);
            free(port_data);
            free(quantized);
        }
        
        // Compute field root from segments
        result->field_roots[field] = gkval_merkle_from_leaves(segments, num_segments);
        result->seg_leaves[field] = segments;
        result->num_seg_leaves[field] = num_segments;
    }
    
    // Compute case root from field roots
    result->case_root = gkval_merkle_from_hashes(result->field_roots, session->manifest.num_fields);
    
    // Check final invariants
    if (!gkval_check_invariants(&result->invariants, &session->manifest.tolerances)) {
        result->success = false;
    }

    // Capture grid info and turning points from reference solver for ledger
    double lam0 = 0.0, dlam = 0.0;
    int tpmi = -1, tpri = -1;
    (void)gkval_geokerr_get_grid_info; // silence unused if not linked
    (void)gkval_geokerr_get_tp_indices;
    // These functions are declared in adapters header and implemented for GeoKerr solver
    gkval_geokerr_get_grid_info(session->geokerr_solver, &lam0, &dlam);
    gkval_geokerr_get_tp_indices(session->geokerr_solver, &tpmi, &tpri);
    // Store invariants and bit patterns
    result->invariants.tpmi = (uint64_t)(tpmi);
    result->invariants.tpri = (uint64_t)(tpri);
    // Bit-cast doubles
    uint64_t lam0_bits = 0, dlam_bits = 0;
    memcpy(&lam0_bits, &lam0, sizeof(uint64_t));
    memcpy(&dlam_bits, &dlam, sizeof(uint64_t));
    result->lam0_bits = lam0_bits;
    result->dlam_bits = dlam_bits;

    // Append ledger line
    char json_line[2048];
    if (gkval_case_result_to_json(result, session->manifest.run_id, session->manifest.grid.S, json_line, sizeof(json_line)) == GKVAL_SUCCESS) {
        gkval_append_ledger_line(session->run_dir, json_line);
    }
    
    // Cleanup solvers
    session->geokerr_solver->cleanup_case(session->geokerr_solver);
    session->port_solver->cleanup_case(session->port_solver);
    
    session->processed_cases++;
    if (!result->success) {
        session->failed_cases++;
    }
    
    return GKVAL_SUCCESS;
}

// =============================================================================
// FAIL WINDOW MANAGEMENT
// =============================================================================

int gkval_create_fail_window(gkval_session_t* session, uint64_t case_id,
                            gkval_field_t field, uint64_t segment_start,
                            uint64_t fail_index_local,
                            const double* geokerr_data, const double* port_data,
                            uint64_t data_len) {
    
    if (!session || !geokerr_data || !port_data || fail_index_local >= data_len) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    // Calculate window bounds
    uint64_t window_size = GKVAL_DEFAULT_WINDOW_SIZE;
    uint64_t start_local = (fail_index_local >= window_size/2) ? (fail_index_local - window_size/2) : 0;
    uint64_t end_local = start_local + window_size;
    if (end_local > data_len) {
        end_local = data_len;
        start_local = (end_local >= window_size) ? (end_local - window_size) : 0;
    }
    uint64_t start_abs = segment_start + start_local;
    uint64_t end_abs = segment_start + end_local;
    
    // Ensure fail_windows/<case_id> directory exists
    char window_path[1024];
    snprintf(window_path, sizeof(window_path), 
             "%s/fail_windows/%lu/%s_%lu_%lu.bin",
             session->run_dir, case_id, gkval_field_to_string(field), start_abs, end_abs);
    
    // Create case subdirectory
    char case_dir[1024];
    snprintf(case_dir, sizeof(case_dir), "%s/fail_windows/%lu", session->run_dir, case_id);
    if (mkdir(case_dir, 0755) != 0 && errno != EEXIST) {
        return GKVAL_ERROR_IO;
    }
    
    return gkval_save_fail_window(session->run_dir, case_id, field, start_abs, end_abs,
                                 geokerr_data + start_local, port_data + start_local);
}

// =============================================================================
// PROGRESS AND CHECKPOINTING
// =============================================================================

int gkval_update_cursor(gkval_session_t* session, uint64_t case_id,
                       gkval_field_t field, uint64_t start) {
    if (!session) return GKVAL_ERROR_INVALID_ARG;
    
    session->cursor.case_id = case_id;
    session->cursor.field = field;
    session->cursor.start = start;
    
    return gkval_save_cursor(session->run_dir, &session->cursor);
}

void gkval_print_progress(const gkval_session_t* session) {
    if (!session) return;
    
    double elapsed = gkval_get_timestamp() - session->start_time;
    double samples_per_sec = (elapsed > 0) ? (session->processed_samples / elapsed) : 0;
    
    printf("Progress: %lu cases, %lu segments, %lu samples (%.1f samples/sec)\n",
           session->processed_cases, session->processed_segments, 
           session->processed_samples, samples_per_sec);
    printf("Failures: %lu cases failed\n", session->failed_cases);
    printf("Current: case %lu, field %s, start %lu\n",
           session->cursor.case_id, 
           gkval_field_to_string(session->cursor.field),
           session->cursor.start);
}

// =============================================================================
// PLACEHOLDER IMPLEMENTATIONS
// =============================================================================

int gkval_validate_cases(gkval_session_t* session, const gkval_inputs_t* cases,
                        uint64_t num_cases, uint64_t start_case_id) {
    // TODO: Implement batch case validation
    return GKVAL_ERROR_INVALID_ARG;
}

int gkval_resume_validation(gkval_session_t* session) {
    // TODO: Implement resume from cursor
    return GKVAL_ERROR_INVALID_ARG;
}

int gkval_checkpoint_progress(gkval_session_t* session) {
    // TODO: Implement progress checkpointing
    return GKVAL_SUCCESS;
}

int gkval_case_result_to_json(const gkval_case_result_t* result,
                             const char* run_id, uint64_t grid_S,
                             char* json_line, size_t max_len) {
    if (!result || !json_line || !run_id) return GKVAL_ERROR_INVALID_ARG;

    // Helper to hex-encode hash
    char field_roots_hex[GKVAL_MAX_FIELDS][GKVAL_HASH_SIZE * 2 + 1];
    for (size_t f = 0; f < GKVAL_MAX_FIELDS; f++) {
        for (int i = 0; i < GKVAL_HASH_SIZE; i++) {
            sprintf(&field_roots_hex[f][i*2], "%02x", result->field_roots[f].bytes[i] & 0xff);
        }
        field_roots_hex[f][GKVAL_HASH_SIZE*2] = '\0';
    }
    char case_root_hex[GKVAL_HASH_SIZE*2 + 1];
    for (int i = 0; i < GKVAL_HASH_SIZE; i++) {
        sprintf(&case_root_hex[i*2], "%02x", result->case_root.bytes[i] & 0xff);
    }
    case_root_hex[GKVAL_HASH_SIZE*2] = '\0';

    // Serialize minimal required fields
    int tpmi_print = (int)result->invariants.tpmi;
    int tpri_print = (int)result->invariants.tpri;
    int written = snprintf(json_line, max_len,
        "{\"run_id\":\"%s\",\"case_id\":%lu,\"success\":%s,"
        "\"grid\":{\"S\":%lu,\"lam0_bits\":%llu,\"dlam_bits\":%llu},"
        "\"tpmi\":%d,\"tpri\":%d,"
        "\"case_root\":\"%s\"}"
        , run_id, result->case_id, result->success ? "true" : "false",
        (unsigned long)grid_S,
        (unsigned long long)result->lam0_bits, (unsigned long long)result->dlam_bits,
        tpmi_print, tpri_print,
        case_root_hex);
    if (written <= 0 || (size_t)written >= max_len) return GKVAL_ERROR_IO;
    return GKVAL_SUCCESS;
}

int gkval_manifest_to_json(const gkval_manifest_t* manifest,
                          char* json_str, size_t max_len) {
    if (!manifest || !json_str || max_len == 0) return GKVAL_ERROR_INVALID_ARG;
    // Serialize a compact JSON representation
    int written = snprintf(json_str, max_len,
        "{\n"
        "  \"version\": \"%s\",\n"
        "  \"run_id\": \"%s\",\n"
        "  \"geokerr_commit\": \"%s\",\n"
        "  \"port_commit\": \"%s\",\n"
        "  \"grid\": { \"type\": \"%s\", \"S\": %lu, \"delta\": %.15e },\n"
        "  \"tolerances\": { \"atol\": %.15e, \"rtol\": %.15e, \"max_ulps\": %lu },\n"
        "  \"quant_lsb\": %u,\n"
        "  \"fields\": [\"%s\", \"%s\", \"%s\", \"%s\", \"%s\"]\n"
        "}\n",
        GKVAL_VERSION,
        manifest->run_id,
        manifest->geokerr_commit,
        manifest->port_commit,
        manifest->grid.type,
        (unsigned long)manifest->grid.S,
        manifest->grid.delta,
        manifest->tolerances.atol,
        manifest->tolerances.rtol,
        (unsigned long)manifest->tolerances.max_ulps,
        manifest->quant_lsb,
        manifest->fields[0], manifest->fields[1], manifest->fields[2], manifest->fields[3], manifest->fields[4]
    );
    return (written > 0 && (size_t)written < max_len) ? GKVAL_SUCCESS : GKVAL_ERROR_IO;
}

int gkval_json_to_manifest(const char* json_str, gkval_manifest_t* manifest) {
    if (!json_str || !manifest) return GKVAL_ERROR_INVALID_ARG;
    memset(manifest, 0, sizeof(*manifest));
    // Very lightweight parser for our own format
    // Extract string fields
    const char* keys_s[] = {"run_id", "geokerr_commit", "port_commit", "type"};
    char* targets_s[] = {manifest->run_id, manifest->geokerr_commit, manifest->port_commit, manifest->grid.type};
    size_t sizes_s[] = {sizeof(manifest->run_id), sizeof(manifest->geokerr_commit), sizeof(manifest->port_commit), sizeof(manifest->grid.type)};
    for (int i = 0; i < 4; i++) {
        const char* k = keys_s[i];
        const char* p = strstr(json_str, k);
        if (p) {
            const char* q = strchr(p, '"');
            if (q) q = strchr(q + 1, '"');
            if (q) {
                const char* r = strchr(q + 1, '"');
                if (r && (size_t)(r - (q + 1)) < sizes_s[i]) {
                    strncpy(targets_s[i], q + 1, r - (q + 1));
                }
            }
        }
    }
    // Extract numeric fields using sscanf on substrings
    unsigned long S = 0, max_ulps = 0; double delta = 0, atol = 0, rtol = 0; unsigned int quant_lsb = 0;
    {
        const char* p = strstr(json_str, "\"S\""); if (p) sscanf(p, "\"S\"%*[^0-9]%lu", &S);
        p = strstr(json_str, "\"delta\""); if (p) sscanf(p, "\"delta\"%*[^0-9.-]%lf", &delta);
        p = strstr(json_str, "\"atol\""); if (p) sscanf(p, "\"atol\"%*[^0-9.-]%lf", &atol);
        p = strstr(json_str, "\"rtol\""); if (p) sscanf(p, "\"rtol\"%*[^0-9.-]%lf", &rtol);
        p = strstr(json_str, "\"max_ulps\""); if (p) sscanf(p, "\"max_ulps\"%*[^0-9]%lu", &max_ulps);
        p = strstr(json_str, "\"quant_lsb\""); if (p) sscanf(p, "\"quant_lsb\"%*[^0-9]%u", &quant_lsb);
    }
    manifest->grid.S = S;
    manifest->grid.delta = delta;
    manifest->tolerances.atol = atol;
    manifest->tolerances.rtol = rtol;
    manifest->tolerances.max_ulps = max_ulps;
    manifest->quant_lsb = (uint8_t)quant_lsb;
    // Default fields (5)
    const char* default_fields[] = {"u","mu","t","phi","affine"};
    manifest->num_fields = 5;
    for (size_t i = 0; i < manifest->num_fields; i++) {
        strncpy(manifest->fields[i], default_fields[i], GKVAL_MAX_FIELD_NAME - 1);
    }
    return GKVAL_SUCCESS;
}
