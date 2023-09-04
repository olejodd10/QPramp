#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "algs.h"
#include "iterable_set.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

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
static double invh[P_DIM][P_DIM];
static double w[C_DIM];
static double g[C_DIM][P_DIM];
static double s[C_DIM][N_DIM];
static double f[P_DIM][N_DIM];

static double ft[N_DIM][P_DIM];
static double neg_w[C_DIM];
static double neg_s[C_DIM][N_DIM];
static double neg_invh[P_DIM][P_DIM]; // Only used once for initial setup
static double neg_invh_f[M_DIM][N_DIM]; // Actually P_DIM x N_DIM, but we only need the M_DIM first rows
static double neg_invh_gt[P_DIM][C_DIM];
static double neg_g_invh[C_DIM][P_DIM];
static double neg_g_invh_gt[C_DIM][C_DIM];
    
static double y[C_DIM];
static double v[C_DIM];
static double invq[C_DIM][C_DIM];
static double x[SIMULATION_TIMESTEPS+1][N_DIM];
static double u[SIMULATION_TIMESTEPS][M_DIM];
static double t[SIMULATION_TIMESTEPS];

static uint8_t setarr1[C_DIM]; 
static ssize_t setarr2[C_DIM]; 
static ssize_t setarr3[C_DIM]; 
static iterable_set_t a_set = {
    .capacity = C_DIM,
    .elements = setarr1,
    .next = setarr2,
    .prev = setarr3,
};

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

    // Other initialization
    timing_reset();
    negate_vector(C_DIM, w, neg_w);
    negate_matrix(C_DIM, N_DIM, s, neg_s);
    transpose(P_DIM, N_DIM, f, ft);
    negate_matrix(P_DIM, P_DIM, invh, neg_invh);
    matrix_product(M_DIM, P_DIM, N_DIM, neg_invh, ft, neg_invh_f);
    matrix_product(P_DIM, P_DIM, C_DIM, neg_invh, g, neg_invh_gt);
    matrix_product(C_DIM, P_DIM, P_DIM, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(C_DIM, P_DIM, C_DIM, g, neg_g_invh, neg_g_invh_gt);
    matrix_product(C_DIM, P_DIM, M_DIM, g, neg_invh, neg_g_invh); // Make sure memory layout is correct for later use
    set_init(&a_set);
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < TEST_CASES; ++i) {
        // Initial state
        memcpy(x[0], x0[i], sizeof(double)*N_DIM);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();
            algorithm2(C_DIM, N_DIM, M_DIM, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_g_invh, x[j], invq, &a_set, y, v, u[j]);
            simulate(N_DIM, M_DIM, a, x[j], b, u[j], x[j+1]); 
            t[j] = (double)timing_elapsed();
            test_case_time += t[j];
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
        if (!matrix_eq(SIMULATION_TIMESTEPS+1, N_DIM, x, x_ref)) {
            printf("WARNING: Verification failed for x in test case %d\n", i);
        }
        fclose(f);

        sprintf(path, U_REF_DIR "/uout%d.csv", i);
        f = fopen(path, "r");
        if (parse_matrix_csv(path, SIMULATION_TIMESTEPS, M_DIM, u_ref)) {
            printf("Error while parsing reference matrix u_ref.\n"); 
        }
        if (!matrix_eq(SIMULATION_TIMESTEPS, M_DIM, u, u_ref)) {
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

    return 0;
}
