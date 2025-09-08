/*
 * GEOKERR VALIDATION CORE LIBRARY IMPLEMENTATION
 * 
 * Implements cryptographic validation with ULP-based quantization,
 * BLAKE3 hashing, and Merkle tree construction.
 */

#define _GNU_SOURCE
#include "gkval_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#ifdef __SSE2__
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

// Simple BLAKE3 implementation (for now using placeholder)
// In production, use the official BLAKE3 C library
static void blake3_hash(const void* data, size_t len, uint8_t* output) {
    // Placeholder: use a simple hash for now
    // TODO: Replace with actual BLAKE3 implementation
    uint64_t hash = 0x1234567890ABCDEF;
    const uint8_t* bytes = (const uint8_t*)data;
    
    for (size_t i = 0; i < len; i++) {
        hash ^= bytes[i];
        hash *= 0x100000001B3ULL;
    }
    
    // Fill output with derived hash
    for (int i = 0; i < GKVAL_HASH_SIZE; i++) {
        output[i] = (uint8_t)(hash >> (8 * (i % 8)));
    }
}

// =============================================================================
// INITIALIZATION
// =============================================================================

int gkval_init(void) {
    // Initialize any global state if needed
    return GKVAL_SUCCESS;
}

int gkval_validate_fp_env(bool fix) {
#if defined(__SSE2__)
    // MXCSR layout: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#techs=SSE
    // Bits of interest:
    //  - Rounding control (RC): bits 13-14 should be 00b (round to nearest)
    //  - Flush-to-zero (FTZ): bit 15 must be 0 (disabled)
    //  - Denormals-are-zero (DAZ): bit 6 must be 0 (disabled)
    unsigned int mxcsr = _mm_getcsr();
    unsigned int desired_mask = 0;
    unsigned int rc_bits = (mxcsr >> 13) & 0x3;
    int ok = (rc_bits == 0) && ((mxcsr & (1u << 15)) == 0) && ((mxcsr & (1u << 6)) == 0);
    if (ok) return GKVAL_SUCCESS;
    if (!fix) return GKVAL_ERROR_IO;

    // Clear RC to 00, clear FTZ and DAZ
    mxcsr &= ~(3u << 13); // RC=00
    mxcsr &= ~(1u << 15); // FTZ=0
    mxcsr &= ~(1u << 6);  // DAZ=0
    _mm_setcsr(mxcsr);

    // Re-read to confirm
    unsigned int mxcsr2 = _mm_getcsr();
    rc_bits = (mxcsr2 >> 13) & 0x3;
    ok = (rc_bits == 0) && ((mxcsr2 & (1u << 15)) == 0) && ((mxcsr2 & (1u << 6)) == 0);
    return ok ? GKVAL_SUCCESS : GKVAL_ERROR_IO;
#else
    // Non-x86/SSE targets: nothing to validate
    (void)fix;
    return GKVAL_SUCCESS;
#endif
}

// =============================================================================
// HASH FUNCTIONS WITH DOMAIN SEPARATION
// =============================================================================

gkval_hash_t gkval_hash_leaf(const char* run_id, uint64_t case_id,
                            gkval_field_t field, uint64_t start, uint64_t end,
                            const uint8_t* quantized_data, size_t data_len) {
    gkval_hash_t result;
    
    // Create domain-separated input: prefix("leaf") || run_id || case_id || field || start || end || data
    size_t prefix_len = strlen("leaf");
    size_t run_id_len = strlen(run_id);
    size_t total_len = prefix_len + run_id_len + sizeof(case_id) + sizeof(field) + 
                      sizeof(start) + sizeof(end) + data_len;
    
    uint8_t* input = malloc(total_len);
    if (!input) {
        memset(result.bytes, 0, GKVAL_HASH_SIZE);
        return result;
    }
    
    size_t offset = 0;
    
    // Add prefix
    memcpy(input + offset, "leaf", prefix_len);
    offset += prefix_len;
    
    // Add run_id
    memcpy(input + offset, run_id, run_id_len);
    offset += run_id_len;
    
    // Add case_id (little-endian)
    uint64_t case_id_le = case_id;  // Assume host is little-endian
    memcpy(input + offset, &case_id_le, sizeof(case_id_le));
    offset += sizeof(case_id_le);
    
    // Add field
    uint32_t field_le = (uint32_t)field;
    memcpy(input + offset, &field_le, sizeof(field_le));
    offset += sizeof(field_le);
    
    // Add start and end
    memcpy(input + offset, &start, sizeof(start));
    offset += sizeof(start);
    memcpy(input + offset, &end, sizeof(end));
    offset += sizeof(end);
    
    // Add quantized data
    memcpy(input + offset, quantized_data, data_len);
    
    // Compute hash
    blake3_hash(input, total_len, result.bytes);
    
    free(input);
    return result;
}

gkval_hash_t gkval_hash_node(const gkval_hash_t* left, const gkval_hash_t* right) {
    gkval_hash_t result;
    
    // Create input: prefix("node") || left_root || right_root
    const char* prefix = "node";
    size_t prefix_len = strlen(prefix);
    size_t total_len = prefix_len + 2 * GKVAL_HASH_SIZE;
    
    uint8_t input[prefix_len + 2 * GKVAL_HASH_SIZE];
    
    memcpy(input, prefix, prefix_len);
    memcpy(input + prefix_len, left->bytes, GKVAL_HASH_SIZE);
    memcpy(input + prefix_len + GKVAL_HASH_SIZE, right->bytes, GKVAL_HASH_SIZE);
    
    blake3_hash(input, total_len, result.bytes);
    return result;
}

// =============================================================================
// ULP-BASED QUANTIZATION
// =============================================================================

int gkval_quantize_double(double x, uint8_t quant_lsb, uint64_t* quantized) {
    if (!quantized || quant_lsb > 52) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    // Reinterpret as uint64
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    
    // Check for NaN/Inf
    uint64_t exp = (bits >> 52) & 0x7FF;
    if (exp == 0x7FF) {
        return GKVAL_ERROR_QUANTIZE;  // NaN/Inf not allowed
    }
    
    if (quant_lsb == 0) {
        *quantized = bits;
        return GKVAL_SUCCESS;
    }
    
    // Extract components
    uint64_t sign = bits & 0x8000000000000000ULL;
    uint64_t frac = bits & 0x000FFFFFFFFFFFFFULL;
    
    // Round to nearest even at bit quant_lsb
    uint64_t mask = (1ULL << quant_lsb) - 1;
    uint64_t half = 1ULL << (quant_lsb - 1);
    uint64_t tail = frac & mask;
    uint64_t frac_q = frac & ~mask;
    
    // Round to nearest, ties to even
    if (tail > half || (tail == half && ((frac_q >> quant_lsb) & 1))) {
        frac_q += (1ULL << quant_lsb);
        
        // Handle carry into exponent
        if (frac_q >= (1ULL << 52)) {
            frac_q = 0;
            exp += 1;
            if (exp >= 0x7FF) {
                return GKVAL_ERROR_QUANTIZE;  // Overflow to Inf
            }
        }
    }
    
    // Reassemble
    *quantized = sign | (exp << 52) | (frac_q & 0x000FFFFFFFFFFFFFULL);
    return GKVAL_SUCCESS;
}

int gkval_quantize_block(const double* data, size_t n, uint8_t quant_lsb,
                        uint8_t* output, size_t* output_len) {
    if (!data || !output || !output_len || *output_len < n * sizeof(uint64_t)) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    uint64_t* quantized = (uint64_t*)output;
    
    for (size_t i = 0; i < n; i++) {
        int ret = gkval_quantize_double(data[i], quant_lsb, &quantized[i]);
        if (ret != GKVAL_SUCCESS) {
            return ret;
        }
    }
    
    *output_len = n * sizeof(uint64_t);
    return GKVAL_SUCCESS;
}

// =============================================================================
// ULP DISTANCE CALCULATION
// =============================================================================

// Map IEEE-754 double to monotonically increasing 64-bit integer space.
// Negative numbers are ordered before positive, with bitwise mapping that
// preserves ordering so that integer difference equals ULP distance.
static inline uint64_t gkval_ordered_uint64(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof(u));
    // If sign bit is set (negative), flip all bits; else flip only the sign bit.
    // This yields a total order compatible with IEEE-754 and two's-complement.
    return (u & 0x8000000000000000ULL) ? ~u : (u ^ 0x8000000000000000ULL);
}

uint64_t gkval_ulp_distance(double a, double b) {
    if (isnan(a) || isnan(b)) return UINT64_MAX;
    if (isinf(a) || isinf(b)) {
        // If both are the same infinity, distance 0; else treat as very large
        return (a == b) ? 0 : UINT64_MAX;
    }
    if (a == b) return 0;

    uint64_t oa = gkval_ordered_uint64(a);
    uint64_t ob = gkval_ordered_uint64(b);
    return (oa > ob) ? (oa - ob) : (ob - oa);
}

// =============================================================================
// ACCEPTANCE PREDICATE
// =============================================================================

bool gkval_accept_sample(double a, double b, const gkval_tolerances_t* tol) {
    if (!tol) return false;
    
    // Check absolute/relative tolerance
    double abs_diff = fabs(a - b);
    double max_val = fmax(fabs(a), fabs(b));
    bool abs_rel_ok = (abs_diff <= tol->atol + tol->rtol * max_val);

    // ULP check
    uint64_t ulp_dist = gkval_ulp_distance(a, b);
    bool ulp_ok = (ulp_dist <= tol->max_ulps);

    // Debug printf can be noisy; keep commented unless troubleshooting
    // printf("ACCEPT_SAMPLE: a=%.10f, b=%.10f, abs_diff=%.2e, max=%.6f, thr=%.2e, abs_rel_ok=%d, ulp=%lu, ulp_ok=%d\n",
    //        a, b, abs_diff, max_val, tol->atol + tol->rtol * max_val, abs_rel_ok, (unsigned long)ulp_dist, ulp_ok);

    return abs_rel_ok && ulp_ok;
}

int gkval_accept_block(const double* a, const double* b, size_t n,
                      const gkval_tolerances_t* tol, size_t* first_fail) {
    if (!a || !b || !tol) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    for (size_t i = 0; i < n; i++) {
        if (!gkval_accept_sample(a[i], b[i], tol)) {
            if (first_fail) *first_fail = i;
            return GKVAL_ERROR_QUANTIZE;  // Reuse for "mismatch"
        }
    }
    
    return GKVAL_SUCCESS;
}

// =============================================================================
// ANGLE UNWRAPPING
// =============================================================================

double gkval_unwrap_angle(double angle_prev, double angle_curr) {
    double diff = angle_curr - angle_prev;
    if (diff > M_PI) {
        return angle_curr - 2.0 * M_PI;
    } else if (diff < -M_PI) {
        return angle_curr + 2.0 * M_PI;
    }
    return angle_curr;
}

int gkval_unwrap_angles(double* angles, size_t n) {
    if (!angles || n == 0) return GKVAL_ERROR_INVALID_ARG;
    
    for (size_t i = 1; i < n; i++) {
        angles[i] = gkval_unwrap_angle(angles[i-1], angles[i]);
    }
    
    return GKVAL_SUCCESS;
}

// =============================================================================
// MERKLE TREE OPERATIONS
// =============================================================================

gkval_hash_t gkval_merkle_from_leaves(const gkval_segment_leaf_t* leaves, size_t n) {
    gkval_hash_t zero_hash = {{0}};
    
    if (!leaves || n == 0) return zero_hash;
    if (n == 1) return leaves[0].root;
    
    // Create array of hashes from leaves
    gkval_hash_t* hashes = malloc(n * sizeof(gkval_hash_t));
    if (!hashes) return zero_hash;
    
    for (size_t i = 0; i < n; i++) {
        hashes[i] = leaves[i].root;
    }
    
    gkval_hash_t result = gkval_merkle_from_hashes(hashes, n);
    free(hashes);
    return result;
}

gkval_hash_t gkval_merkle_from_hashes(const gkval_hash_t* hashes, size_t n) {
    gkval_hash_t zero_hash = {{0}};
    
    if (!hashes || n == 0) return zero_hash;
    if (n == 1) return hashes[0];
    
    // Build Merkle tree bottom-up
    gkval_hash_t* current_level = malloc(n * sizeof(gkval_hash_t));
    if (!current_level) return zero_hash;
    
    memcpy(current_level, hashes, n * sizeof(gkval_hash_t));
    size_t current_n = n;
    
    while (current_n > 1) {
        size_t next_n = (current_n + 1) / 2;
        
        for (size_t i = 0; i < next_n; i++) {
            if (2 * i + 1 < current_n) {
                // Pair exists
                current_level[i] = gkval_hash_node(&current_level[2*i], &current_level[2*i+1]);
            } else {
                // Odd node, promote directly
                current_level[i] = current_level[2*i];
            }
        }
        current_n = next_n;
    }
    
    gkval_hash_t result = current_level[0];
    free(current_level);
    return result;
}

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

const char* gkval_field_to_string(gkval_field_t field) {
    switch (field) {
        case GKVAL_FIELD_U: return "u";
        case GKVAL_FIELD_MU: return "mu"; 
        case GKVAL_FIELD_T: return "t";
        case GKVAL_FIELD_PHI: return "phi";
        case GKVAL_FIELD_AFFINE: return "affine";
        default: return "unknown";
    }
}

gkval_field_t gkval_string_to_field(const char* str) {
    if (!str) return GKVAL_FIELD_U;
    if (strcmp(str, "u") == 0) return GKVAL_FIELD_U;
    if (strcmp(str, "mu") == 0) return GKVAL_FIELD_MU;
    if (strcmp(str, "t") == 0) return GKVAL_FIELD_T;
    if (strcmp(str, "phi") == 0) return GKVAL_FIELD_PHI;
    if (strcmp(str, "affine") == 0) return GKVAL_FIELD_AFFINE;
    return GKVAL_FIELD_U;
}

void gkval_hash_to_hex(const gkval_hash_t* hash, char* hex_str) {
    if (!hash || !hex_str) return;
    
    for (int i = 0; i < GKVAL_HASH_SIZE; i++) {
        sprintf(hex_str + 2*i, "%02x", hash->bytes[i]);
    }
    hex_str[2 * GKVAL_HASH_SIZE] = '\0';
}

int gkval_hex_to_hash(const char* hex_str, gkval_hash_t* hash) {
    if (!hex_str || !hash || strlen(hex_str) != 2 * GKVAL_HASH_SIZE) {
        return GKVAL_ERROR_INVALID_ARG;
    }
    
    for (int i = 0; i < GKVAL_HASH_SIZE; i++) {
        unsigned int byte;
        if (sscanf(hex_str + 2*i, "%02x", &byte) != 1) {
            return GKVAL_ERROR_INVALID_ARG;
        }
        hash->bytes[i] = (uint8_t)byte;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_create_run_dir(const char* run_dir) {
    if (!run_dir) return GKVAL_ERROR_INVALID_ARG;
    // Ensure parent directory ("runs") exists for relative paths like "runs/ID"
    // Best-effort: create "runs" if run_dir starts with "runs/"
    if (strncmp(run_dir, "runs/", 5) == 0) {
        if (mkdir("runs", 0755) != 0 && errno != EEXIST) {
            return GKVAL_ERROR_IO;
        }
    }
    if (mkdir(run_dir, 0755) != 0 && errno != EEXIST) {
        return GKVAL_ERROR_IO;
    }
    
    return GKVAL_SUCCESS;
}

int gkval_create_fail_window_dir(const char* run_dir) {
    if (!run_dir) return GKVAL_ERROR_INVALID_ARG;
    
    char fail_dir[512];
    snprintf(fail_dir, sizeof(fail_dir), "%s/fail_windows", run_dir);
    
    if (mkdir(fail_dir, 0755) != 0 && errno != EEXIST) {
        return GKVAL_ERROR_IO;
    }
    
    return GKVAL_SUCCESS;
}

const char* gkval_strerror(int error_code) {
    switch (error_code) {
        case GKVAL_SUCCESS: return "Success";
        case GKVAL_ERROR_INVALID_ARG: return "Invalid argument";
        case GKVAL_ERROR_IO: return "I/O error";
        case GKVAL_ERROR_HASH: return "Hash error";
        case GKVAL_ERROR_QUANTIZE: return "Quantization error";
        case GKVAL_ERROR_MERKLE: return "Merkle tree error";
        case GKVAL_ERROR_JSON: return "JSON error";
        case GKVAL_ERROR_MEMORY: return "Memory error";
        default: return "Unknown error";
    }
}