#include <stdio.h>
#include <stdint.h>

#include "algs.h"
#include "vector.h"
#include "matrix.h"
#include "parse.h"
#include "lti.h"

#define SIMULATION_TIMESTEPS 10

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

static double x[N];
static double a[N][N];
static double b[N][M];

static double invh[P][P];
static double w[C];
static double g[C][P];
static double s[C][N];
static double f[P][N];

static double ft[N][P];
static double neg_w[C];
static double neg_s[C][N];
static double neg_invh[P][P]; // Only used once for initial setup
static double neg_invh_f[P][N];
static double neg_invh_gt[P][C];
static double neg_g_invh[C][P];
static double neg_g_invh_gt[C][C];
    
static double y[C];
static double v[C];
static double invq[C][C];
static uint8_t a_set[C]; // Lookup table to represent set for now
static double temp1[C];
static double temp2[C][C];
static double temp3[P];
static double temp4[N];
    
static double u[P];

int main() {
    // Parsing input dimensions
    // The buffers probably take a lot of memory
    FILE* f_a = fopen(A_PATH, "r");
    FILE* f_b = fopen(B_PATH, "r");
    FILE* f_x0 = fopen(X0_PATH, "r");
    FILE* f_invh = fopen(INVH_PATH, "r");
    FILE* f_w = fopen(W_PATH, "r");
    FILE* f_g = fopen(G_PATH, "r");
    FILE* f_s = fopen(S_PATH, "r");
    FILE* f_f = fopen(F_PATH, "r");

    // Parsing input vectors and input matrices
    // TODO: Assert their width and height first
    parse_matrix_csv(f_a, N, N, a);
    parse_matrix_csv(f_b, N, M, b);
    parse_vector_csv(f_x0, N, x);
    parse_matrix_csv(f_invh, P, P, invh);
    parse_vector_csv(f_w, C, w);
    parse_matrix_csv(f_g, C, P, g);
    parse_matrix_csv(f_s, C, N, s);
    parse_matrix_csv(f_f, P, N, f);
    fclose(f_a);
    fclose(f_b);
    fclose(f_x0);
    fclose(f_invh);
    fclose(f_w);
    fclose(f_g);
    fclose(f_s);
    fclose(f_f);

    // Other initialization
    negate_vector(C, w, neg_w);
    negate_matrix(C, N, s, neg_s);
    transpose(P, N, f, ft);
    negate_matrix(P, P, invh, neg_invh);
    matrix_product(P, P, N, neg_invh, ft, neg_invh_f);
    matrix_product(P, P, C, neg_invh, g, neg_invh_gt);
    matrix_product(C, P, P, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(C, P, C, g, neg_g_invh, neg_g_invh_gt);

    // Simulation
    for (uint16_t i = 0; i < SIMULATION_TIMESTEPS; ++i) {
        algorithm2(C, N, P, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_invh_gt, x, invq, a_set, y, v, temp1, temp2, temp3, u);
        simulate(N, M, a, x, b, u, temp4, x); // Note that the second dimension of u is larger than M
    }
    printf("Simulation finished with the following state vector:\n");
    print_vector(N, x);

    return 0;
}
