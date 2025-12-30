#include "gkval_validator.h"
#include "gkval_core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Stub adapter accessors to satisfy linker when pulling in gkval_validator.o.
int gkval_geokerr_get_grid_info(const gkval_solver_t* solver, double* lam0, double* dlam) {
    (void)solver;
    if (lam0) *lam0 = 0.0;
    if (dlam) *dlam = 0.0;
    return GKVAL_SUCCESS;
}

int gkval_geokerr_get_tp_indices(const gkval_solver_t* solver, int* tpmi, int* tpri) {
    (void)solver;
    if (tpmi) *tpmi = 0;
    if (tpri) *tpri = 0;
    return GKVAL_SUCCESS;
}

static int check_fail_window_contents(const double* geokerr_out, const double* port_out,
                                      size_t n, size_t start_local) {
    for (size_t i = 0; i < n; ++i) {
        double expect_geokerr = 1000.0 + (double)(start_local + i);
        double expect_port = 2000.0 + (double)(start_local + i);
        if (geokerr_out[i] != expect_geokerr || port_out[i] != expect_port) {
            fprintf(stderr, "Mismatch at i=%zu: got %.1f/%.1f, expected %.1f/%.1f\n",
                    i, geokerr_out[i], port_out[i], expect_geokerr, expect_port);
            return 1;
        }
    }
    return 0;
}

int main(void) {
    const char* run_dir = "runs/FAILWIN_TEST";
    const uint64_t case_id = 1;
    const gkval_field_t field = GKVAL_FIELD_PHI;
    const uint64_t segment_start = 1024;
    const uint64_t data_len = 1024;
    const uint64_t fail_index_local = 100;

    if (gkval_create_run_dir(run_dir) != GKVAL_SUCCESS) {
        fprintf(stderr, "Failed to create run dir: %s\n", run_dir);
        return 1;
    }
    if (gkval_create_fail_window_dir(run_dir) != GKVAL_SUCCESS) {
        fprintf(stderr, "Failed to create fail window dir: %s\n", run_dir);
        return 1;
    }

    gkval_session_t session;
    memset(&session, 0, sizeof(session));
    strncpy(session.run_dir, run_dir, sizeof(session.run_dir) - 1);

    double* geokerr_data = (double*)malloc(data_len * sizeof(double));
    double* port_data = (double*)malloc(data_len * sizeof(double));
    if (!geokerr_data || !port_data) {
        fprintf(stderr, "Allocation failed\n");
        free(geokerr_data);
        free(port_data);
        return 1;
    }

    for (uint64_t i = 0; i < data_len; ++i) {
        geokerr_data[i] = 1000.0 + (double)i;
        port_data[i] = 2000.0 + (double)i;
    }

    if (gkval_create_fail_window(&session, case_id, field, segment_start, fail_index_local,
                                 geokerr_data, port_data, data_len) != GKVAL_SUCCESS) {
        fprintf(stderr, "gkval_create_fail_window failed\n");
        free(geokerr_data);
        free(port_data);
        return 1;
    }

    uint64_t window_size = GKVAL_DEFAULT_WINDOW_SIZE;
    uint64_t start_local = (fail_index_local >= window_size / 2)
                           ? (fail_index_local - window_size / 2)
                           : 0;
    uint64_t end_local = start_local + window_size;
    if (end_local > data_len) {
        end_local = data_len;
        start_local = (end_local >= window_size) ? (end_local - window_size) : 0;
    }
    uint64_t start_abs = segment_start + start_local;
    uint64_t end_abs = segment_start + end_local;
    size_t out_len = (size_t)(end_local - start_local);

    double* geokerr_out = (double*)malloc(out_len * sizeof(double));
    double* port_out = (double*)malloc(out_len * sizeof(double));
    if (!geokerr_out || !port_out) {
        fprintf(stderr, "Output allocation failed\n");
        free(geokerr_data);
        free(port_data);
        free(geokerr_out);
        free(port_out);
        return 1;
    }

    if (gkval_load_fail_window(run_dir, case_id, field, start_abs, end_abs,
                               geokerr_out, port_out) != GKVAL_SUCCESS) {
        fprintf(stderr, "gkval_load_fail_window failed\n");
        free(geokerr_data);
        free(port_data);
        free(geokerr_out);
        free(port_out);
        return 1;
    }

    if (check_fail_window_contents(geokerr_out, port_out, out_len, (size_t)start_local) != 0) {
        free(geokerr_data);
        free(port_data);
        free(geokerr_out);
        free(port_out);
        return 1;
    }

    printf("PASS: fail window stored with global indices [%llu, %llu)\n",
           (unsigned long long)start_abs, (unsigned long long)end_abs);

    free(geokerr_data);
    free(port_data);
    free(geokerr_out);
    free(port_out);
    return 0;
}
