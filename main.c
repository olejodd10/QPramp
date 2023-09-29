#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

#include "osqp.h"
#include "cs.h"
#include "auxil.h"
#include "osqp_init.h"

#define SIMULATION_TIMESTEPS 100

#define INPUT_DIR "../examples/osqp/n20"
#define OUTPUT_DIR INPUT_DIR "/out"
#define REFERENCE_DIR INPUT_DIR "/reference"

#define A_PATH INPUT_DIR "/a.csv"
#define B_PATH INPUT_DIR "/b.csv"
#define X0_PATH INPUT_DIR "/x0.csv"

#define OSQP_P_PATH INPUT_DIR "/osqp_p.csv"
#define OSQP_Q_PATH INPUT_DIR "/osqp_q.csv"
#define OSQP_L_PATH INPUT_DIR "/osqp_l.csv"
#define OSQP_A_PATH INPUT_DIR "/osqp_a.csv"
#define OSQP_U_PATH INPUT_DIR "/osqp_u.csv"

#define X_DIR OUTPUT_DIR
#define U_DIR OUTPUT_DIR
#define T_DIR OUTPUT_DIR

#ifdef REFERENCE_DIR
#define X_REF_DIR REFERENCE_DIR
#define U_REF_DIR REFERENCE_DIR
#define EQ_EPS 1e-8
#endif

static void print_csc(csc *mat) {
    printf("%lld %lld %lld %lld \n", mat->nzmax, mat->m, mat->n, mat->nz);
}

int main() {
    timing_print_precision();
    timing_reset();
	size_t initial_conditions = csv_parse_matrix_height(X0_PATH);
	size_t n_dim = csv_parse_matrix_width(A_PATH);
	size_t m_dim = csv_parse_matrix_width(B_PATH);
	size_t c_dim = csv_parse_matrix_height(OSQP_A_PATH);
	size_t p_dim  = csv_parse_matrix_width(OSQP_A_PATH);
    size_t horizon = (p_dim - n_dim)/(n_dim + m_dim);
    printf("Input dimension parsing time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    double *a = (double*)malloc(n_dim*n_dim*sizeof(double));
    double *b = (double*)malloc(n_dim*m_dim*sizeof(double));
    double *x0 = (double*)malloc(initial_conditions*n_dim*sizeof(double));

	double *osqp_p = (double*)malloc(p_dim*p_dim*sizeof(double));
	double *osqp_q = (double*)malloc(p_dim*sizeof(double));
	double *osqp_l = (double*)malloc(c_dim*sizeof(double));
	double *osqp_a = (double*)malloc(c_dim*p_dim*sizeof(double));
	double *osqp_u = (double*)malloc(c_dim*sizeof(double));

    double *x = (double*)malloc((SIMULATION_TIMESTEPS+1)*n_dim*sizeof(double));
    double *u = (double*)malloc(SIMULATION_TIMESTEPS*m_dim*sizeof(double));
	double *t = (double*)malloc(SIMULATION_TIMESTEPS*sizeof(double));

#ifdef REFERENCE_DIR
    double *x_ref = (double*)malloc((SIMULATION_TIMESTEPS+1)*n_dim*sizeof(double));
    double *u_ref = (double*)malloc(SIMULATION_TIMESTEPS*m_dim*sizeof(double));
#endif
    printf("Allocation time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    if (csv_parse_matrix(A_PATH, n_dim, n_dim, a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    if (csv_parse_matrix(B_PATH, n_dim, m_dim, b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    if (csv_parse_matrix(X0_PATH, initial_conditions, n_dim, x0)) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }

    if (csv_parse_matrix(OSQP_P_PATH, p_dim, p_dim, osqp_p)) { 
        printf("Error while parsing input matrix osqp_p.\n");
        return 1;
    }
    if (csv_parse_vector(OSQP_Q_PATH, p_dim, osqp_q)) { 
        printf("Error while parsing input vector osqp_q.\n");
        return 1;
    }
    if (csv_parse_vector(OSQP_L_PATH, c_dim, osqp_l)) { 
        printf("Error while parsing input vector osqp_l.\n");
        return 1;
    }
    if (csv_parse_matrix(OSQP_A_PATH, c_dim, p_dim, osqp_a)) { 
        printf("Error while parsing input matrix osqp_a.\n");
        return 1;
    }
    if (csv_parse_vector(OSQP_U_PATH, c_dim, osqp_u)) { 
        printf("Error while parsing input vector osqp_u.\n");
        return 1;
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    const size_t nnz_a = osqp_init_num_nonzero(c_dim, p_dim, osqp_a);
    const size_t nnz_p = osqp_init_num_nonzero_symmetric(p_dim, osqp_p);
    c_float p_x[nnz_p];
    c_int p_i[nnz_p];
    c_int p_p[p_dim+1];
    c_float a_x[nnz_a];
    c_int a_i[nnz_a];
    c_int a_p[p_dim+1];
	memset(p_x, 0, sizeof(c_float)*nnz_p);
	memset(p_i, 0, sizeof(c_int)*nnz_p);
	memset(p_p, 0, sizeof(c_int)*(p_dim+1));
	memset(a_x, 0, sizeof(c_float)*nnz_a);
	memset(a_i, 0, sizeof(c_int)*nnz_a);
	memset(a_p, 0, sizeof(c_int)*(p_dim+1));

    // Other initialization
    timing_reset();
    // printf("Dette er nnz: %d %d\n", osqp_init_num_nonzero(p_dim, p_dim, h), osqp_init_num_nonzero(c_dim, p_dim, g));

    // Populate constant matrices and vectors
    osqp_init_populate_csc(c_dim, p_dim, nnz_a, osqp_a, a_x, a_i, a_p);
    csc *osqp_a_csc = csc_matrix(c_dim, p_dim, nnz_a, a_x, a_i, a_p);
    osqp_init_populate_upper_triangular_csc(p_dim, nnz_p, osqp_p, p_x, p_i, p_p);
    csc *osqp_p_csc = csc_matrix(p_dim, p_dim, nnz_p, p_x, p_i, p_p);
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    /* Solver, settings, matrices */
    OSQPWorkspace   *workspace;
    OSQPData data = (OSQPData) {
        .n = p_dim,
        .m = c_dim,
        .P = osqp_p_csc,
        .A = osqp_a_csc,
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
    }

    /* Setup workspace */
    c_int exitflag = osqp_setup(&workspace, &data, settings);
    // printf("Exitflag after setup: %lld\n", exitflag);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < initial_conditions; ++i) {
        // Initial state
        memcpy(&x[0], &x0[i*n_dim], sizeof(double)*n_dim);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();

            /* Solve problem */
            exitflag = osqp_solve(workspace);
            // printf("Exitflag after solve: %lld\n", exitflag);
            double *u_temp = (double*)(&workspace->solution->x[n_dim*(horizon+1)]);

            // Apply input
            simulate(n_dim, m_dim, a, &x[j*n_dim], b, u_temp, &x[(j+1)*n_dim]); 

            // Update constraints
            for (uint16_t k = 0; k < n_dim; ++k) {
                osqp_l[k] = -x[(j+1)*n_dim + k];
                osqp_u[k] = -x[(j+1)*n_dim + k];
            }
            osqp_update_bounds(workspace, osqp_l, osqp_u);

            t[j] = (double)timing_elapsed();
            test_case_time += t[j];

            memcpy(&u[j*m_dim], u_temp, m_dim*sizeof(double)); // Save for output

            // printf("Solution u: ");
            // for (uint16_t k = 0; k < m_dim; ++k) {
            //     printf("%.2e ", u[j*m_dim+k]);
            // }
            // printf("\n");

            // printf("New x: ");
            // for (uint16_t k = 0; k < n_dim; ++k) {
            //     printf("%.2e ", x[j*n_dim+k]);
            // }
            // printf("\n");
        }
        total_time += test_case_time;

        // Test case result printing
        // printf("Total simulation time for %d iterations of test case %d: %.0f us\n", SIMULATION_TIMESTEPS, i, test_case_time/1000);
        // print_vector(n_dim, x[SIMULATION_TIMESTEPS]);

        char path[80];
        #ifdef REFERENCE_DIR
        // Verify test case results
        sprintf(path, X_REF_DIR "/xout%d.csv", i);
        FILE* f = fopen(path, "r");
        if (csv_parse_matrix(path, SIMULATION_TIMESTEPS+1, n_dim, x_ref)) {
            printf("Error while parsing reference matrix x_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS+1, n_dim, x, x_ref, EQ_EPS)) {
            printf("WARNING: Verification failed for x in test case %d\n", i);
        }
        fclose(f);

        sprintf(path, U_REF_DIR "/uout%d.csv", i);
        f = fopen(path, "r");
        if (csv_parse_matrix(path, SIMULATION_TIMESTEPS, m_dim, u_ref)) {
            printf("Error while parsing reference matrix u_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS, m_dim, u, u_ref, EQ_EPS)) {
            printf("WARNING: Verification failed for u in test case %d\n", i);
        }
        fclose(f);
        #endif

        // Save test case results
        sprintf(path, X_DIR "/xout%d.csv", i);
        if (csv_save_matrix(path, SIMULATION_TIMESTEPS+1, n_dim, x) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(path, U_DIR "/uout%d.csv", i);
        if (csv_save_matrix(path, SIMULATION_TIMESTEPS, m_dim, u) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(path, T_DIR "/tout%d.csv", i);
        if (csv_save_vector(path, SIMULATION_TIMESTEPS, t) < 0) {
            printf("Error while saving t.\n");
        }
    }

    // Summary
    printf("Average time per test case running %d timesteps: %.0f us\n", SIMULATION_TIMESTEPS, total_time/1000/initial_conditions);

    /* Cleanup */
    osqp_cleanup(workspace);
    if (settings) free(settings);

    return 0;
}
