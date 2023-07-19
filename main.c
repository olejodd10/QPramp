#include <stdio.h>
#include <stdint.h>

#include "algs.h"
#include "vector.h"
#include "matrix.h"
#include "parse.h"
#include "lti.h"

#define SIMULATION_TIMESTEPS 10

#define EXAMPLE_PATH "../examples/example1"

#define A_PATH EXAMPLE_PATH "/a.csv"
#define B_PATH EXAMPLE_PATH "/b.csv"
#define X0_PATH EXAMPLE_PATH "/x0.csv"
#define INVH_PATH EXAMPLE_PATH "/invh.csv"
#define W_PATH EXAMPLE_PATH "/w.csv"
#define G_PATH EXAMPLE_PATH "/g.csv"
#define S_PATH EXAMPLE_PATH "/s.csv"
#define F_PATH EXAMPLE_PATH "/f.csv"

int main() {
    // Parsing input dimensions
    FILE* f_a = fopen(A_PATH, "r");
    FILE* f_b = fopen(B_PATH, "r");
    FILE* f_x0 = fopen(X0_PATH, "r");
    FILE* f_invh = fopen(INVH_PATH, "r");
    FILE* f_w = fopen(W_PATH, "r");
    FILE* f_g = fopen(G_PATH, "r");
    FILE* f_s = fopen(S_PATH, "r");
    FILE* f_f = fopen(F_PATH, "r");

    size_t n = parse_matrix_csv_width(f_s);
    size_t c = parse_matrix_csv_height(f_s);
    size_t p = parse_matrix_csv_width(f_g);
    size_t m = parse_matrix_csv_width(f_b);

    // Allocation
    double x[n];
    double a[n][n];
    double b[n][m];

    double invh[p][p];
    double w[c];
    double g[c][p];
    double s[c][n];
    double f[p][n];
    
    double ft[n][p];
    double neg_w[c];
    double neg_s[c][n];
    double neg_invh[p][p]; // Only used once for initial setup
    double neg_invh_f[p][n];
    double neg_invh_gt[p][c];
    double neg_g_invh[c][p];
    double neg_g_invh_gt[c][c];
    
    double y[c];
    double v[c];
    double invq[c][c];
    uint8_t a_set[c]; // Lookup table to represent set for now
    double temp1[c];
    double temp2[c][c];
    double temp3[p];
    double temp4[n];
    
    double u[p];

    // Parsing input vectors and input matrices
    parse_matrix_csv(f_a, n, n, a);
    parse_matrix_csv(f_b, n, m, b);
    parse_vector_csv(f_x0, n, x);
    parse_matrix_csv(f_invh, p, p, invh);
    parse_vector_csv(f_w, c, w);
    parse_matrix_csv(f_g, c, p, g);
    parse_matrix_csv(f_s, c, n, s);
    parse_matrix_csv(f_f, p, n, f);
    fclose(f_a);
    fclose(f_b);
    fclose(f_x0);
    fclose(f_invh);
    fclose(f_w);
    fclose(f_g);
    fclose(f_s);
    fclose(f_f);

    // Other initialization
    negate_vector(c, w, neg_w);
    negate_matrix(c, n, s, neg_s);
    transpose(p, n, f, ft);
    negate_matrix(p, p, invh, neg_invh);
    matrix_product(p, p, n, neg_invh, ft, neg_invh_f);
    matrix_product(p, p, c, neg_invh, g, neg_invh_gt);
    matrix_product(c, p, p, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(c, p, c, g, neg_g_invh, neg_g_invh_gt);

    // Simulation
    for (uint16_t i = 0; i < SIMULATION_TIMESTEPS; ++i) {
        algorithm2(c, n, p, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_invh_gt, x, invq, a_set, y, v, temp1, temp2, temp3, u);
        simulate(n, m, a, x, b, u, temp4, x); // Note that the second dimension of u is larger than m
    }
    printf("Simulation finished with the following state vector:\n");
    print_vector(n, x);

    return 0;
}
