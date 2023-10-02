#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "qp_ramp.h"
#include "iterable_set.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

#ifndef SIMULATION_TIMESTEPS
#error "SIMULATION_TIMESTEPS not set"
#endif

#ifndef INPUT_DIR
#error "INPUT_DIR not set"
#endif

#define A_PATH INPUT_DIR "/a.csv"
#define B_PATH INPUT_DIR "/b.csv"
#define X0_PATH INPUT_DIR "/x0.csv"
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
#endif

int main() {
    timing_print_precision();
    timing_reset();
	size_t initial_conditions = csv_parse_matrix_height(X0_PATH);
	size_t n_dim = csv_parse_matrix_width(A_PATH);
	size_t m_dim = csv_parse_matrix_width(B_PATH);
	size_t c_dim = csv_parse_matrix_height(G_PATH);
	size_t p_dim  = csv_parse_matrix_width(G_PATH);
    printf("Input dimension parsing time: %ld us\n", timing_elapsed()/1000);

    timing_reset();
    double *a = (double*)malloc(n_dim*n_dim*sizeof(double));
    double *b = (double*)malloc(n_dim*m_dim*sizeof(double));
    double *x0 = (double*)malloc(initial_conditions*n_dim*sizeof(double));
    double *invh = (double*)malloc(p_dim*p_dim*sizeof(double));
	double *w = (double*)malloc(c_dim*sizeof(double));
    double *g = (double*)malloc(c_dim*p_dim*sizeof(double));
    double *s = (double*)malloc(c_dim*n_dim*sizeof(double));
    double *f = (double*)malloc(p_dim*n_dim*sizeof(double));

    double *ft = (double*)malloc(n_dim*p_dim*sizeof(double));
	double *neg_w = (double*)malloc(c_dim*sizeof(double));
    double *neg_s = (double*)malloc(c_dim*n_dim*sizeof(double));
    double *neg_invh = (double*)malloc(p_dim*p_dim*sizeof(double)); // Only used once for initial setup
    double *neg_invh_f = (double*)malloc(m_dim*n_dim*sizeof(double)); // Actually p_dim x n_dim, but we only need the m_dim first rows
    double *neg_invh_gt = (double*)malloc(p_dim*c_dim*sizeof(double));
    double *neg_g_invh = (double*)malloc(c_dim*p_dim*sizeof(double));
    double *neg_g_invh_gt = (double*)malloc(c_dim*c_dim*sizeof(double));
        
	double *y = (double*)malloc(c_dim*sizeof(double));
	double *v = (double*)malloc(c_dim*sizeof(double));
    double *invq = (double*)malloc(c_dim*c_dim*sizeof(double));
    double *x = (double*)malloc((SIMULATION_TIMESTEPS+1)*n_dim*sizeof(double));
    double *u = (double*)malloc(SIMULATION_TIMESTEPS*m_dim*sizeof(double));
	double *t = (double*)malloc(SIMULATION_TIMESTEPS*sizeof(double));

	uint8_t *setarr1 = (uint8_t*)malloc(c_dim*sizeof(uint8_t)); 
	ssize_t *setarr2 = (ssize_t*)malloc(c_dim*sizeof(ssize_t)); 
	ssize_t *setarr3 = (ssize_t*)malloc(c_dim*sizeof(ssize_t)); 
    iterable_set_t a_set = {
        .capacity = c_dim,
        .elements = setarr1,
        .next = setarr2,
        .prev = setarr3,
    };

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
    if (csv_parse_matrix(INVH_PATH, p_dim, p_dim, invh)) { 
        printf("Error while parsing input matrix invh.\n"); 
        return 1; 
    }
    if (csv_parse_vector(W_PATH, c_dim, w)) { 
        printf("Error while parsing input vector w.\n"); 
        return 1; 
    }
    if (csv_parse_matrix(G_PATH, c_dim, p_dim, g)) { 
        printf("Error while parsing input matrix g.\n"); 
        return 1; 
    }
    if (csv_parse_matrix(S_PATH, c_dim, n_dim, s)) { 
        printf("Error while parsing input matrix s.\n"); 
        return 1; 
    }
    if (csv_parse_matrix(F_PATH, p_dim, n_dim, f)) { 
        printf("Error while parsing input matrix f.\n"); 
        return 1; 
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    // Other initialization
    timing_reset();
    negate_vector(c_dim, w, neg_w);
    negate_matrix(c_dim, n_dim, s, neg_s);
    transpose(p_dim, n_dim, f, ft);
    negate_matrix(p_dim, p_dim, invh, neg_invh);
    matrix_product(m_dim, p_dim, n_dim, neg_invh, ft, neg_invh_f);
    matrix_product(p_dim, p_dim, c_dim, neg_invh, g, neg_invh_gt);
    matrix_product(c_dim, p_dim, p_dim, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(c_dim, p_dim, c_dim, g, neg_g_invh, neg_g_invh_gt);
    matrix_product(c_dim, p_dim, m_dim, g, neg_invh, neg_g_invh); // Make sure memory layout is correct for later use
    set_init(&a_set);
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < initial_conditions; ++i) {
        // Initial state
        memcpy(&x[0], &x0[i*n_dim], sizeof(double)*n_dim);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();
            qp_ramp_solve_mpc(c_dim, n_dim, m_dim, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_g_invh, &x[j*n_dim], invq, &a_set, y, v, &u[j*m_dim]);
            simulate(n_dim, m_dim, a, &x[j*n_dim], b, &u[j*m_dim], &x[(j+1)*n_dim]); 
            t[j] = (double)timing_elapsed();
            test_case_time += t[j];
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
        if (!matrix_eq(SIMULATION_TIMESTEPS+1, n_dim, x, x_ref, EPS)) {
            printf("WARNING: Verification failed for x in test case %d\n", i);
        }
        fclose(f);

        sprintf(path, U_REF_DIR "/uout%d.csv", i);
        f = fopen(path, "r");
        if (csv_parse_matrix(path, SIMULATION_TIMESTEPS, m_dim, u_ref)) {
            printf("Error while parsing reference matrix u_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS, m_dim, u, u_ref, EPS)) {
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

    return 0;
}
