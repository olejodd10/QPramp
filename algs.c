#include "algs.h"

static uint8_t compatible(uint16_t c, uint8_t a_set[c], double y[c]) {
    for (int i = 0; i < c; ++i) {
        if ((a_set[i] && y[i] < 0.0) || (!a_set[i] && y[i] >= 0.0)) {
            return 0;
        }
    }
    return 1;
}

static int most_negative_index(uint16_t c, uint8_t a_set[c], double y[c]) {
    double min = y[0];
    int index = 0;
    for (int i = 1; i < c; ++i) {
        if (a_set[i] && y[i] < min) {
            min = y[i];
            index = i;
        }    
    }
    if (!a_set[index] || y[index] > 0.0) {
        return -1;
    }
    return index;
}

static int most_positive_index(uint16_t c, uint8_t a_set[c], double y[c]) {
    double max = y[0];
    int index = 0;
    for (int i = 1; i < c; ++i) {
        if (!a_set[i] && y[i] > max) {
            max = y[i];
            index = i;
        }    
    }
    if (a_set[index] || y[index] < 0.0) {
        return -1;
    }
    return index;
}

// ******** Remember that neg_g_invh_gt must be volumn major for this implementation to work! Or maybe it is symmetric?
static void compute_v(uint16_t c, double invq[c][c], uint8_t a_contains_i, uint16_t i, double neg_g_invh_gt[c][c], double v[c]) {
    vector_copy(c, neg_g_invh_gt[i], v); // Ensure column major!
    if (!a_contains_i) {
        v[i] += 1.0;
    }
    matrix_vector_product(c, c, invq, v, v);
}

// Terminal results are written to and available from a_set, invq and y
// v can be anything - call ite temp?
static void algorithm1(uint16_t c, double invq[c][c], uint8_t a_set[c], double neg_g_invh_gt[c][c], double v[c], double temp1[c], double temp2[c][c], double y[c]) {
    int index;
    double q0;
    while (!compatible(c, a_set, y)) {
        // Tror faktisk ikke rustfunksjonene jeg skrev gjør noe som er raskere?
        index = most_negative_index(c, a_set, y);
        if (index >= 0) {
            q0 = 1.0;
            compute_v(c, invq, 1, index, neg_g_invh_gt, v);
            a_set[index] = 0;
        } else {
            index = most_positive_index(c, a_set, y);
            if (index < 0) {
                printf("Don't think this should happen");
            }
            q0 = -1.0;
            compute_v(c, invq, 0, index, neg_g_invh_gt, v);
            a_set[index] = 1;
        }
        
        // Now that we have q0 and v, compute invq and y
        scale_vector(c, v, -1.0/(q0+v[index]), v);
        vector_copy(c, invq[index], temp1);
        outer_product(c, c, v, temp1, temp2); // **Sparsity
        matrix_sum(c, c, invq, temp2, invq); // **Sparsity
        scale_vector(c, v, y[index], v);
        vector_sum(c, y, v, y);
    }
}
//
// Changes y
static void ramp(uint16_t c, double y[c]) {
    for (int i = 0; i < c; ++i) {
        if (y[i] < 0.0) {
            y[i] = 0.0;
        }
    }
}

static void compute_u(uint16_t p, uint16_t n, uint16_t c, double neg_invh_f[p][n], double x[n], double neg_invh_gt[p][c], double y[c], double temp[p], double u[p]) {
    matrix_vector_product(p, n, neg_invh_f, x, u);
    ramp(c, y);
    matrix_vector_product(p, c, neg_invh_gt, y, temp); // **Sparsity
    vector_sum(p, u, temp, u);
}

// This only represents one iteration of algorithm2, i.e. algorithm2 without the lines with "while"
// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
// Note that y is modified
void algorithm2(uint16_t c, uint16_t n, uint16_t p, double neg_g_invh_gt[c][c], double neg_s[c][n], double neg_w[c], double neg_invh_f[p][n], double neg_invh_gt[p][c], double x[n], double invq[c][c], uint8_t a_set[c], double y[c], double v[c], double temp1[c], double temp2[c][c], double temp3[p], double u[p]) {
    eye(c, invq);
    for (int i = 0; i < c; ++i) { //Later: call clear set function
        a_set[i] = 0;
    }
    matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, invq, a_set, neg_g_invh_gt, v, temp1, temp2, y);
    compute_u(p, n, c, neg_invh_f, x, neg_invh_gt, y, temp3, u);
}
