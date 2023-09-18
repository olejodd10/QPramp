#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

#include "osqp.h"
#include "cs.h"
#include "auxil.h"
#include "osqp_init.h"

#ifndef TEST_CASES
#error "TEST_CASES not set"
#endif

#ifndef SIMULATION_TIMESTEPS
#error "SIMULATION_TIMESTEPS not set"
#endif

#ifndef N_DIM
#error "N_DIM not set"
#endif

#ifndef M_DIM
#error "M_DIM not set"
#endif

#ifndef HORIZON
#error "HORIZON not set"
#endif

#ifndef C_DIM
#error "C_DIM not set"
#endif

#define P_DIM HORIZON*M_DIM

#ifndef INPUT_DIR
#error "INPUT_DIR not set"
#endif

#define A_PATH INPUT_DIR "/a.csv"
#define B_PATH INPUT_DIR "/b.csv"
#define X0_PATH INPUT_DIR "/x0.csv"
#define H_PATH INPUT_DIR "/h.csv"
#define INVH_PATH INPUT_DIR "/invh.csv"
#define W_PATH INPUT_DIR "/w.csv"
#define G_PATH INPUT_DIR "/g.csv"
#define S_PATH INPUT_DIR "/s.csv"
#define F_PATH INPUT_DIR "/f.csv"

#define X_DIR OUTPUT_DIR
#define U_DIR OUTPUT_DIR
#define T_DIR OUTPUT_DIR

#ifdef REFERENCE_DIR
#define X_REF_DIR REFERENCE_DIR
#define U_REF_DIR REFERENCE_DIR
static double x_ref[SIMULATION_TIMESTEPS+1][N_DIM];
static double u_ref[SIMULATION_TIMESTEPS][M_DIM];
#endif

static double a[N_DIM][N_DIM];
static double b[N_DIM][M_DIM];
static double x0[TEST_CASES][N_DIM];
static double h[P_DIM][P_DIM];
static double invh[P_DIM][P_DIM];
static double w[C_DIM];
static double g[C_DIM][P_DIM];
static double s[C_DIM][N_DIM];
static double f[P_DIM][N_DIM];

static double ft[N_DIM][P_DIM];
static double invh_f[P_DIM][N_DIM];

static c_int P_p[P_DIM+1];
static c_int A_p[P_DIM+1];
static c_float osqp_l[C_DIM];
static c_float osqp_u[C_DIM];
static c_float osqp_q[P_DIM];

static double x[SIMULATION_TIMESTEPS+1][N_DIM];
static double u[SIMULATION_TIMESTEPS][M_DIM];
static double t[SIMULATION_TIMESTEPS];

static void print_csc(csc *mat) {
    printf("%lld %lld %lld %lld \n", mat->nzmax, mat->m, mat->n, mat->nz);
}

int main() {
    timing_print_precision();
    timing_reset();
    if (parse_matrix_csv(A_PATH, N_DIM, N_DIM, a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    if (parse_matrix_csv(B_PATH, N_DIM, M_DIM, b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    if (parse_matrix_csv(X0_PATH, TEST_CASES, N_DIM, x0)) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }
    if (parse_matrix_csv(H_PATH, P_DIM, P_DIM, h)) { 
        printf("Error while parsing input matrix h.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(INVH_PATH, P_DIM, P_DIM, invh)) { 
        printf("Error while parsing input matrix invh.\n"); 
        return 1; 
    }
    if (parse_vector_csv(W_PATH, C_DIM, w)) { 
        printf("Error while parsing input vector w.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(G_PATH, C_DIM, P_DIM, g)) { 
        printf("Error while parsing input matrix g.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(S_PATH, C_DIM, N_DIM, s)) { 
        printf("Error while parsing input matrix s.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(F_PATH, P_DIM, N_DIM, f)) { 
        printf("Error while parsing input matrix f.\n"); 
        return 1; 
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    const size_t nnz_g = osqp_init_num_nonzero(C_DIM, P_DIM, g);
    const size_t nnz_h = osqp_init_num_nonzero_symmetric(P_DIM, h);
    c_float P_x[nnz_h];
    c_int P_i[nnz_h];
    c_float A_x[nnz_g];
    c_int A_i[nnz_g];
	memset(P_x, 0, sizeof(c_float)*nnz_h);
	memset(P_i, 0, sizeof(c_int)*nnz_h);
	memset(P_p, 0, sizeof(c_int)*(P_DIM+1));
	memset(A_x, 0, sizeof(c_float)*nnz_g);
	memset(A_i, 0, sizeof(c_int)*nnz_g);
	memset(A_p, 0, sizeof(c_int)*(P_DIM+1));
	memset(osqp_l, 0, sizeof(c_float)*C_DIM);
	memset(osqp_u, 0, sizeof(c_float)*C_DIM);
	memset(osqp_q, 0, sizeof(c_float)*P_DIM);

    // Other initialization
    timing_reset();
    transpose(P_DIM, N_DIM, f, ft);
    matrix_product(P_DIM, P_DIM, N_DIM, invh, ft, invh_f);

    // printf("Dette er nnz: %d %d\n", osqp_init_num_nonzero(P_DIM, P_DIM, h), osqp_init_num_nonzero(C_DIM, P_DIM, g));

    // Populate constant matrices and vectors
    osqp_init_populate_csc(C_DIM, P_DIM, nnz_g, g, A_x, A_i, A_p);
    csc *osqp_a = csc_matrix(C_DIM, P_DIM, nnz_g, A_x, A_i, A_p);
    osqp_init_populate_upper_triangular_csc(P_DIM, nnz_h, h, P_x, P_i, P_p);
    csc *osqp_p = csc_matrix(P_DIM, P_DIM, nnz_h, P_x, P_i, P_p);
    for (size_t i = 0; i < C_DIM; ++i) {
        osqp_l[i] = -OSQP_INFTY;
    }
    memset(osqp_u, 0, sizeof(c_float)*P_DIM); // Overwritten later
    memset(osqp_q, 0, sizeof(c_float)*P_DIM);
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    /* Solver, settings, matrices */
    OSQPWorkspace   *workspace;
    OSQPData data = (OSQPData) {
        .n = P_DIM,
        .m = C_DIM,
        .P = osqp_p,
        .A = osqp_a,
        .q = osqp_q,
        .l = osqp_l,
        .u = osqp_u,
    };
    OSQPSettings    *settings;

    /* Set default settings */
    settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    if (settings) {
        osqp_set_default_settings(settings);
        settings->verbose = 0;
        settings->check_termination = 2;
        settings->warm_start = 1;
        settings->polish = 0;
    }

    /* Setup workspace */
    c_int exitflag = osqp_setup(&workspace, &data, settings);
    // printf("Exitflag after setup: %lld\n", exitflag);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < TEST_CASES; ++i) {
        // Initial state
        memcpy(x[0], x0[i], sizeof(double)*N_DIM);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();

            // Update upper constraints
            matrix_vector_product(C_DIM, N_DIM, s, x[j], osqp_u);
            vector_sum(C_DIM, osqp_u, w, osqp_u);
            exitflag = osqp_update_upper_bound(workspace, osqp_u);
            // printf("Exitflag after upper bound update: %lld\n", exitflag);

            /* Solve problem */
            exitflag = osqp_solve(workspace);
            // printf("Exitflag after solve: %lld\n", exitflag);
            double *z = (double*)workspace->solution->x;
            for (uint16_t k = 0; k < M_DIM; ++k) {
                u[j][k] = z[k] - inner_product(N_DIM, invh_f[k], x[j]);
            }

            // Apply input
            simulate(N_DIM, M_DIM, a, x[j], b, u[j], x[j+1]); 

            t[j] = (double)timing_elapsed();
            test_case_time += t[j];

            // printf("Solution z: ");
            // for (uint16_t k = 0; k < P_DIM; ++k) {
            //     printf("%.2e ", z[k]);
            // }
            // printf("\n");

            // printf("Solution u: ");
            // for (uint16_t k = 0; k < M_DIM; ++k) {
            //     printf("%.2e ", u[j][k]);
            // }
            // printf("\n");

            // printf("New x: ");
            // for (uint16_t k = 0; k < N_DIM; ++k) {
            //     printf("%.2e ", x[j][k]);
            // }
            // printf("\n");
        }
        total_time += test_case_time;

        // Test case result printing
        // printf("Total simulation time for %d iterations of test case %d: %.0f us\n", SIMULATION_TIMESTEPS, i, test_case_time/1000);
        // print_vector(N_DIM, x[SIMULATION_TIMESTEPS]);

        char path[80];
        #ifdef REFERENCE_DIR
        // Verify test case results
        sprintf(path, X_REF_DIR "/xout%d.csv", i);
        FILE* f = fopen(path, "r");
        if (parse_matrix_csv(path, SIMULATION_TIMESTEPS+1, N_DIM, x_ref)) {
            printf("Error while parsing reference matrix x_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS+1, N_DIM, x, x_ref, EPS)) {
            printf("WARNING: Verification failed for x in test case %d\n", i);
        }
        fclose(f);

        sprintf(path, U_REF_DIR "/uout%d.csv", i);
        f = fopen(path, "r");
        if (parse_matrix_csv(path, SIMULATION_TIMESTEPS, M_DIM, u_ref)) {
            printf("Error while parsing reference matrix u_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS, M_DIM, u, u_ref, EPS)) {
            printf("WARNING: Verification failed for u in test case %d\n", i);
        }
        fclose(f);
        #endif

        // Save test case results
        sprintf(path, X_DIR "/xout%d.csv", i);
        if (save_matrix_csv(path, SIMULATION_TIMESTEPS+1, N_DIM, x) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(path, U_DIR "/uout%d.csv", i);
        if (save_matrix_csv(path, SIMULATION_TIMESTEPS, M_DIM, u) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(path, T_DIR "/tout%d.csv", i);
        if (save_vector_csv(path, SIMULATION_TIMESTEPS, t) < 0) {
            printf("Error while saving t.\n");
        }
    }

    // Summary
    printf("Average time per test case running %d timesteps: %.0f us\n", SIMULATION_TIMESTEPS, total_time/1000/TEST_CASES);

    /* Cleanup */
    osqp_cleanup(workspace);
    if (settings) free(settings);

    return 0;
}
