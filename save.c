#include <stdio.h>
#include <stdint.h>

#include "vector.h"
#include "matrix.h"
#include "csv.h"
#include "timing.h"

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

#define X0_PATH INPUT_DIR "/x0.csv"
#define INVH_PATH INPUT_DIR "/invh.csv"
#define W_PATH INPUT_DIR "/w.csv"
#define G_PATH INPUT_DIR "/g.csv"
#define S_PATH INPUT_DIR "/s.csv"
#define F_PATH INPUT_DIR "/f.csv"

#define NEG_W_PATH INPUT_DIR "/neg_w.csv"
#define NEG_S_PATH INPUT_DIR "/neg_s.csv"
#define NEG_INVH_F_PATH INPUT_DIR "/neg_invh_f.csv"
#define NEG_G_INVH_PATH INPUT_DIR "/neg_g_invh.csv"
#define NEG_G_INVH_GT_PATH INPUT_DIR "/neg_g_invh_gt.csv"

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

int main() {
    tick();
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

    tick();
    negate_vector(C, w, neg_w);
    negate_matrix(C, N, s, neg_s);
    transpose(P, N, f, ft);
    negate_matrix(P, P, invh, neg_invh);
    matrix_product(M, P, N, neg_invh, ft, neg_invh_f);
    matrix_product(P, P, C, neg_invh, g, neg_invh_gt);
    matrix_product(C, P, P, g, neg_invh, neg_g_invh); // Exploiting the fact that invh is symmetric
    matrix_product(C, P, C, g, neg_g_invh, neg_g_invh_gt);
    printf("Initialization time: %d ms\n", tock());

    // Save matrices
    tick();
	save_vector_csv(NEG_W_PATH, C, neg_w);
	save_matrix_csv(NEG_S_PATH, C, N, neg_s);
	save_matrix_csv(NEG_INVH_F_PATH, P, N, neg_invh_f);
	save_matrix_csv(NEG_G_INVH_PATH, C, P, neg_g_invh);
    save_matrix_csv(NEG_G_INVH_GT_PATH, C, C, neg_g_invh_gt);
    printf("Save time: %d ms\n", tock());

    return 0;
}
