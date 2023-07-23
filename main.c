#include <stdio.h>
#include <stdint.h>

#include "algs.h"
#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "lti.h"
#include "timing.h"

#define SIMULATION_TIMESTEPS 100

#define EXAMPLE 3

#if EXAMPLE == 1
    #define N 2
    #define M 1
    #define HORIZON 10
    #define C 66 
#elif EXAMPLE == 2
    #define N 4
    #define M 2
    #define HORIZON 30
    #define C 316 
#elif EXAMPLE == 3
    #define N 25
    #define M 4
    #define HORIZON 30
    #define C 1560 
#endif
#define P HORIZON*M

#define SWAP(a,b) { typeof(a) SWAP = a; a = b; b = SWAP; }

#define STR(x) STR2(x)
#define STR2(x) #x
#define EXAMPLE_PATH "../examples/example" STR(EXAMPLE)

#define A_PATH EXAMPLE_PATH "/a.csv"
#define B_PATH EXAMPLE_PATH "/b.csv"
#define X0_PATH EXAMPLE_PATH "/x0.csv"
#define INVH_PATH EXAMPLE_PATH "/invh.csv"
#define W_PATH EXAMPLE_PATH "/w.csv"
#define G_PATH EXAMPLE_PATH "/g.csv"
#define S_PATH EXAMPLE_PATH "/s.csv"
#define F_PATH EXAMPLE_PATH "/f.csv"

static double a[N][N];
static double b[N][M];
static double x0[N];
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
static uint8_t a_set[C]; // Lookup table to represent set for now
static double x1[N];
    
static double u[M];

int main() {
    tick();
    if (parse_matrix_csv(A_PATH, N, N, a)) { 
        printf("Error while parsing input matrix a.\n"); 
        return 1;
    }
    if (parse_matrix_csv(B_PATH, N, M, b)) { 
        printf("Error while parsing input matrix b.\n"); 
        return 1;
    }
    if (parse_vector_csv(X0_PATH, N, x0)) { 
        printf("Error while parsing input vector x0.\n"); 
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
    printf("Input parsing time: %d ms\n", tock());

    // Other initialization
    tick();
    negate_vector(C, w, neg_w);
    negate_matrix(C, N, s, neg_s);
    transpose(P, N, f, ft);
    negate_matrix(P, P, invh, neg_invh);
    matrix_product(M, P, N, neg_invh, ft, neg_invh_f);
    matrix_product(P, P, C, neg_invh, g, neg_invh_gt);
    matrix_product(C, P, P, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(C, P, C, g, neg_g_invh, neg_g_invh_gt);
    matrix_product(C, P, M, g, neg_invh, neg_g_invh); // Make sure memory layout is correct for later use
    printf("Initialization time: %d ms\n", tock());

    // Simulation
    double* x0_p = x0;
    double* x1_p = x1;
    tick();
    for (uint16_t i = 0; i < SIMULATION_TIMESTEPS; ++i) {
        algorithm2(C, N, M, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_g_invh, x0_p, invq, a_set, y, v, u);
        simulate(N, M, a, x0_p, b, u, x1_p); 
        SWAP(x0_p, x1_p);
    }
    printf("Simulation time for %d iterations: %d ms\n", SIMULATION_TIMESTEPS, tock());
    printf("Simulation finished with the following state vector:\n");
    print_vector(N, x0);

    return 0;
}
