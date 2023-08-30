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

#ifndef N
#error "N not set"
#endif

#ifndef M
#error "M not set"
#endif

#ifndef HORIZON
#error "HORIZON not set"
#endif

#ifndef C
#error "C not set"
#endif

#define P HORIZON*M

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

static double a[N][N];
static double b[N][M];
static double x0[TEST_CASES][N];
static double invh[P][P];
static double w[C];
static double g[C][P];
static double s[C][N];
static double f[P][N];

static double ft[N][P];
static double neg_w[C];
static double neg_s[C][N];
static double neg_invh[P][P]; // Only used once for initial setup
static double neg_invh_f[M][N]; // Actually P x N, but we only need the M first rows
static double neg_invh_gt[P][C];
static double neg_g_invh[C][P];
static double neg_g_invh_gt[C][C];
    
static double y[C];
static double v[C];
static double invq[C][C];
static double x[SIMULATION_TIMESTEPS+1][N];
static double u[SIMULATION_TIMESTEPS][M];
static double t[SIMULATION_TIMESTEPS];

static uint8_t setarr1[C]; 
static ssize_t setarr2[C]; 
static ssize_t setarr3[C]; 
static iterable_set_t a_set = {
    .capacity = C,
    .elements = setarr1,
    .next = setarr2,
    .prev = setarr3,
};

int main() {
    timing_print_precision();
    timing_reset();
    if (parse_matrix_csv(A_PATH, N, N, a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    if (parse_matrix_csv(B_PATH, N, M, b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    if (parse_matrix_csv(X0_PATH, TEST_CASES, N, x0)) { 
        printf("Error while parsing input matrix x0.\n"); 
        return 1;
    }
    if (parse_matrix_csv(INVH_PATH, P, P, invh)) { 
        printf("Error while parsing input matrix invh.\n"); 
        return 1; 
    }
    if (parse_vector_csv(W_PATH, C, w)) { 
        printf("Error while parsing input vector w.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(G_PATH, C, P, g)) { 
        printf("Error while parsing input matrix g.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(S_PATH, C, N, s)) { 
        printf("Error while parsing input matrix s.\n"); 
        return 1; 
    }
    if (parse_matrix_csv(F_PATH, P, N, f)) { 
        printf("Error while parsing input matrix f.\n"); 
        return 1; 
    }
    printf("Input parsing time: %ld us\n", timing_elapsed()/1000);

    // Other initialization
    timing_reset();
    negate_vector(C, w, neg_w);
    negate_matrix(C, N, s, neg_s);
    transpose(P, N, f, ft);
    negate_matrix(P, P, invh, neg_invh);
    matrix_product(M, P, N, neg_invh, ft, neg_invh_f);
    matrix_product(P, P, C, neg_invh, g, neg_invh_gt);
    matrix_product(C, P, P, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(C, P, C, g, neg_g_invh, neg_g_invh_gt);
    matrix_product(C, P, M, g, neg_invh, neg_g_invh); // Make sure memory layout is correct for later use
    set_init(&a_set);
    printf("Initialization time: %ld us\n", timing_elapsed()/1000);

    // Simulation
    double total_time = 0.0;
    for (uint16_t i = 0; i < TEST_CASES; ++i) {
        // Initial state
        memcpy(x[0], x0[i], sizeof(double)*N);

        double test_case_time = 0.0;
        for (uint16_t j = 0; j < SIMULATION_TIMESTEPS; ++j) {
            timing_reset();
            algorithm2(C, N, M, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_g_invh, x[j], invq, &a_set, y, v, u[j]);
            simulate(N, M, a, x[j], b, u[j], x[j+1]); 
            t[j] = (double)timing_elapsed();
            test_case_time += t[j];
        }
        total_time += test_case_time;

        // Test case result printing
        printf("Total simulation time for %d iterations of test case %d: %.0f us\n", SIMULATION_TIMESTEPS, i, test_case_time/1000);
        print_vector(N, x[SIMULATION_TIMESTEPS]);

        // Save test case results
        char str[80];
        sprintf(str, X_DIR "/xout%d.csv", i);
        if (save_matrix_csv(str, SIMULATION_TIMESTEPS+1, N, x) < 0) {
            printf("Error while saving x.\n");
        }
        sprintf(str, U_DIR "/uout%d.csv", i);
        if (save_matrix_csv(str, SIMULATION_TIMESTEPS, M, u) < 0) {
            printf("Error while saving u.\n");
        }
        sprintf(str, T_DIR "/tout%d.csv", i);
        if (save_vector_csv(str, SIMULATION_TIMESTEPS, t) < 0) {
            printf("Error while saving t.\n");
        }
    }

    // Summary
    printf("Average time per test case running %d timesteps: %.0f us\n", SIMULATION_TIMESTEPS, total_time/1000/TEST_CASES);

    return 0;
}
