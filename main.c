#include <stdio.h>
#include <stdint.h>

#include "algs.h"
#include "vector.h"
#include "matrix.h"

#define N 3
#define P 2
#define C 1

static double x[N] = {0.0, -1.0, 4.0};

static double invh[P][P] = {{1.0, 0.0}, {0.0, 1.0}};
static double w[C] = {0.5};
static double g[C][P] = {{0.5, 0.5}};
static double s[C][N] = {{0.9, 0.6, 0.4}};
static double f[P][N] = {{0.5, 0.4, 0.6}, {0.1, 0.7, 0.5}};

// Everything below is set during setup
static double ft[N][P];
static double gt[P][C];
static double neg_w[C];
static double neg_s[C][N];
static double neg_invh[P][P]; // Only used once for initial setup
static double neg_g_invh_gt[C][C];
static double neg_invh_f[P][N];
static double neg_invh_gt[P][C];

static double y[C];
static double v[C];
static double invq[C][C];
static uint8_t a_set[C]; // Lookup table to represent set for now
static double temp1[C];
static double temp2[C][C];
static double temp3[P];

static double u[P];

int main() {
    negate_vector(C, w, neg_w);
    negate_matrix(C, N, s, neg_s);
    transpose(P, N, f, ft);
    transpose(C, P, g, gt);
    negate_matrix(P, P, invh, neg_invh);
    matrix_product(P, P, N, neg_invh, f, neg_invh_f);
    matrix_product(P, P, C, neg_invh, gt, neg_invh_gt);
    matrix_product(C, P, C, g, neg_invh_gt, neg_g_invh_gt);

    // while (1) {
        algorithm2(C, N, P, neg_g_invh_gt, neg_s, neg_w, neg_invh_f, neg_invh_gt, x, invq, a_set, y, v, temp1, temp2, temp3, u);
        printf("Alg 2 returned, this is the resulting u:\n");
        print_vector(P, u);
        // Simulate
    // }
    return 0;
}
