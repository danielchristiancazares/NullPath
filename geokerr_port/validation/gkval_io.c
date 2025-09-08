/*
 * GEOKERR VALIDATION I/O FUNCTIONS
 * 
 * File operations for manifest, cursor, ledger, and fail windows.
 */

#include "gkval_core.h"
#include "gkval_validator.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

// =============================================================================
// MANIFEST OPERATIONS
// =============================================================================

int gkval_save_manifest(const char* run_dir, const gkval_manifest_t* manifest) {
    if (!run_dir || !manifest) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char manifest_path[1024];
    snprintf(manifest_path, sizeof(manifest_path), "%s/manifest.json", run_dir);
    
    char buf[2048];
    if (gkval_manifest_to_json(manifest, buf, sizeof(buf)) != GKVAL_SUCCESS) {
        return GKVAL_ERROR_JSON;
    }
    return gkval_write_file_atomic(manifest_path, buf);
}

int gkval_load_manifest(const char* run_dir, gkval_manifest_t* manifest) {
    if (!run_dir || !manifest) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char manifest_path[1024];
    snprintf(manifest_path, sizeof(manifest_path), "%s/manifest.json", run_dir);
    
    // Read entire file
    FILE* f = fopen(manifest_path, "r");
    if (!f) return GKVAL_ERROR_IO;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz <= 0 || sz > 1<<20) { fclose(f); return GKVAL_ERROR_IO; }
    char* buf = (char*)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return GKVAL_ERROR_MEMORY; }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) { free(buf); fclose(f); return GKVAL_ERROR_IO; }
    buf[sz] = '\0';
    fclose(f);

    int ret = gkval_json_to_manifest(buf, manifest);
    free(buf);
    return ret;
}

// =============================================================================
// CURSOR OPERATIONS
// =============================================================================

int gkval_save_cursor(const char* run_dir, const gkval_cursor_t* cursor) {
    if (!run_dir || !cursor) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char cursor_path[1024];
    snprintf(cursor_path, sizeof(cursor_path), "%s/cursor.json", run_dir);
    
    char content[512];
    snprintf(content, sizeof(content),
             "{\"case_id\": %lu, \"field\": %d, \"start\": %lu}\n",
             cursor->case_id, (int)cursor->field, cursor->start);
    
    return gkval_write_file_atomic(cursor_path, content);
}

int gkval_load_cursor(const char* run_dir, gkval_cursor_t* cursor) {
    if (!run_dir || !cursor) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char cursor_path[1024];
    snprintf(cursor_path, sizeof(cursor_path), "%s/cursor.json", run_dir);
    
    FILE* f = fopen(cursor_path, "r");
    if (!f) return GKVAL_ERROR_IO;
    
    int field_int;
    int ret = fscanf(f, "{\"case_id\": %lu, \"field\": %d, \"start\": %lu}",
                     &cursor->case_id, &field_int, &cursor->start);
    
    fclose(f);
    
    if (ret != 3) return GKVAL_ERROR_IO;
    
    cursor->field = (gkval_field_t)field_int;
    return GKVAL_SUCCESS;
}

// =============================================================================
// LEDGER OPERATIONS
// =============================================================================

int gkval_append_ledger_line(const char* run_dir, const char* json_line) {
    if (!run_dir || !json_line) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char ledger_path[1024];
    snprintf(ledger_path, sizeof(ledger_path), "%s/ledger.ndjson", run_dir);
    
    FILE* f = fopen(ledger_path, "a");
    if (!f) return GKVAL_ERROR_IO;
    
    if (fprintf(f, "%s\n", json_line) < 0) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    if (fflush(f) != 0 || fsync(fileno(f)) != 0) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    fclose(f);
    return GKVAL_SUCCESS;
}

// =============================================================================
// FAIL WINDOW OPERATIONS
// =============================================================================

int gkval_save_fail_window(const char* run_dir, uint64_t case_id,
                          gkval_field_t field, uint64_t start, uint64_t end,
                          const double* geokerr_data, const double* port_data) {
    if (!run_dir || !geokerr_data || !port_data || end <= start) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char window_path[1024];
    snprintf(window_path, sizeof(window_path),
             "%s/fail_windows/%lu/%s_%lu_%lu.bin",
             run_dir, case_id, gkval_field_to_string(field), start, end);
    
    FILE* f = fopen(window_path, "wb");
    if (!f) return GKVAL_ERROR_IO;
    
    // Write header
    uint32_t magic = GKVAL_FAIL_WINDOW_MAGIC;
    uint32_t version = GKVAL_FAIL_WINDOW_VERSION;
    uint16_t field_id = (uint16_t)field;
    uint64_t size = end - start;
    
    if (fwrite(&magic, sizeof(magic), 1, f) != 1 ||
        fwrite(&version, sizeof(version), 1, f) != 1 ||
        fwrite(&field_id, sizeof(field_id), 1, f) != 1 ||
        fwrite(&start, sizeof(start), 1, f) != 1 ||
        fwrite(&end, sizeof(end), 1, f) != 1 ||
        fwrite(&size, sizeof(size), 1, f) != 1) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    // Write data
    if (fwrite(geokerr_data, sizeof(double), size, f) != size ||
        fwrite(port_data, sizeof(double), size, f) != size) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    fclose(f);
    return GKVAL_SUCCESS;
}

int gkval_load_fail_window(const char* run_dir, uint64_t case_id,
                          gkval_field_t field, uint64_t start, uint64_t end,
                          double* geokerr_data, double* port_data) {
    if (!run_dir || !geokerr_data || !port_data || end <= start) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    char window_path[1024];
    snprintf(window_path, sizeof(window_path),
             "%s/fail_windows/%lu/%s_%lu_%lu.bin",
             run_dir, case_id, gkval_field_to_string(field), start, end);
    
    FILE* f = fopen(window_path, "rb");
    if (!f) return GKVAL_ERROR_IO;
    
    // Read and verify header
    uint32_t magic, version;
    uint16_t field_id;
    uint64_t file_start, file_end, file_size;
    
    if (fread(&magic, sizeof(magic), 1, f) != 1 ||
        fread(&version, sizeof(version), 1, f) != 1 ||
        fread(&field_id, sizeof(field_id), 1, f) != 1 ||
        fread(&file_start, sizeof(file_start), 1, f) != 1 ||
        fread(&file_end, sizeof(file_end), 1, f) != 1 ||
        fread(&file_size, sizeof(file_size), 1, f) != 1) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    if (magic != GKVAL_FAIL_WINDOW_MAGIC ||
        version != GKVAL_FAIL_WINDOW_VERSION ||
        field_id != (uint16_t)field ||
        file_start != start || file_end != end) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    uint64_t size = end - start;
    
    // Read data
    if (fread(geokerr_data, sizeof(double), size, f) != size ||
        fread(port_data, sizeof(double), size, f) != size) {
        fclose(f);
        return GKVAL_ERROR_IO;
    }
    
    fclose(f);
    return GKVAL_SUCCESS;
}