#include "algs.h"

#ifndef EPS
#error "EPS not set"
#endif

static ssize_t most_negative_index(size_t c, const uint8_t a_set[c], double y[c]) {
    double min = y[0];
    size_t index = 0;
    for (size_t i = 1; i < c; ++i) {
        if (a_set[i] && y[i] < min) {
            min = y[i];
            index = i;
        }    
    }
    if (!a_set[index] || y[index] > -EPS) {
        return -1;
    }
    return index;
}

static ssize_t most_positive_index(size_t c, const uint8_t a_set[c], double y[c]) {
    double max = y[0];
    size_t index = 0;
    for (size_t i = 1; i < c; ++i) {
        if (y[i] > max && !a_set[i]) {
            max = y[i];
            index = i;
        }    
    }
    if (a_set[index] || y[index] < EPS) {
        return -1;
    }
    return index;
}

static void compute_v(size_t c, const double invq[c][c], const uint8_t a_set[c], size_t index, const double neg_g_invh_gt[c][c], double v[c]) {
    // Compute matrix vector product
    // Sparse part
    for (size_t i = 0; i < c; ++i) {
        v[i] = a_set[i] ? 0.0 : (i == index ? neg_g_invh_gt[index][i] + 1.0 : neg_g_invh_gt[index][i]); // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = 0; i < c; ++i) {
        if (a_set[i]) {
            add_scaled_vector(c, v, invq[i], neg_g_invh_gt[index][i], v);
        }
    }
}

// Terminal results are written to and available from a_set, invq and y
// v can be anything - call it temp?
static void algorithm1(size_t c, double invq[c][c], uint8_t a_set[c], const double neg_g_invh_gt[c][c], double v[c], double y[c]) {
    while (1) {
        double q0 = 1.0;
        ssize_t index = most_negative_index(c, a_set, y);
        if (index >= 0) {
            compute_v(c, invq, a_set, index, neg_g_invh_gt, v);
            a_set[index] = 0;
        } else {
            index = most_positive_index(c, a_set, y);
            if (index < 0) {
                // printf("y and A are compatible\n");
                break;
            }
            compute_v(c, invq, a_set, index, neg_g_invh_gt, v);
            a_set[index] = 1;
            q0 = -1.0;
        }
        
        // Update invq
        scale_vector(c, v, -1.0/(q0+v[index]), v);
        for (size_t i = 0; i < c; ++i) { // Adding the outer product column for column
            if (a_set[i]) {
                if (i == index) { // Just inserted
                    memcpy(invq[i], v, c*sizeof(double));
                    invq[i][i] += 1.0; // Pretend there was a unit vector in the column to start with
                } else {
                    add_scaled_vector(c, invq[i], v, invq[i][index], invq[i]);
                }
            }
        }
        // Update y
        scale_vector(c, v, y[index], v);
        vector_sum(c, y, v, y);
    }
}

static void compute_u(size_t m, size_t n, size_t c, const double neg_invh_f[m][n], const double x[n], const double neg_g_invh[c][m], const double y[c], double u[m]) {
    matrix_vector_product(m, n, neg_invh_f, x, u);
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > EPS) {
            add_scaled_vector(m, u, neg_g_invh[i], y[i], u);
        }
    }
}

// This only represents one iteration of algorithm2, i.e. algorithm2 without the lines with "while"
// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
// Note that y is modified
void algorithm2(size_t c, size_t n, size_t m, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_invh_f[m][n], const double neg_g_invh[c][m], const double x[n], double invq[c][c], uint8_t a_set[c], double y[c], double v[c], double u[m]) {
    memset(a_set, 0, c);
    matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, invq, a_set, neg_g_invh_gt, v, y);
    compute_u(m, n, c, neg_invh_f, x, neg_g_invh, y, u);
}
