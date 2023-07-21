#include "algs.h"

#define EPS 2.2204e-15

static ssize_t most_negative_index(size_t c, uint8_t a_set[c], double y[c]) {
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

static ssize_t most_positive_index(size_t c, uint8_t a_set[c], double y[c]) {
    double max = y[0];
    size_t index = 0;
    for (size_t i = 1; i < c; ++i) {
        if (!a_set[i] && y[i] > max) {
            max = y[i];
            index = i;
        }    
    }
    if (a_set[index] || y[index] < EPS) {
        return -1;
    }
    return index;
}

static void compute_v(size_t c, double invq[c][c], uint8_t a_set[c], size_t index, double neg_g_invh_gt[c][c], double temp[c], double v[c]) {
    vector_copy(c, neg_g_invh_gt[index], temp); // This row is the same as the desired column since the matrix is symmetric
    if (!a_set[index]) {
        temp[index] += 1.0;
    }
    // Compute matrix vector product
    // Sparse part
    for (size_t i = 0; i < c; ++i) {
        v[i] = a_set[i] ? 0.0 : temp[i]; // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = 0; i < c; ++i) {
        if (a_set[i]) {
            add_scaled_vector(c, v, invq[i], temp[i], v);
        }
    }
}

// Terminal results are written to and available from a_set, invq and y
// v can be anything - call it temp?
static void algorithm1(size_t c, double invq[c][c], uint8_t a_set[c], double neg_g_invh_gt[c][c], double v[c], double temp1[c], double y[c]) {
    while (1) {
        double q0 = 1.0;
        ssize_t index = most_negative_index(c, a_set, y);
        if (index >= 0) {
            compute_v(c, invq, a_set, index, neg_g_invh_gt, temp1, v);
            a_set[index] = 0;
        } else {
            index = most_positive_index(c, a_set, y);
            if (index < 0) {
                // printf("y and A are compatible\n");
                break;
            }
            compute_v(c, invq, a_set, index, neg_g_invh_gt, temp1, v);
            a_set[index] = 1;
            q0 = -1.0;
        }
        
        // Update invq
        scale_vector(c, v, -1.0/(q0+v[index]), v);
        for (size_t i = 0; i < c; ++i) { // Adding the outer product column for column
            if (a_set[i]) {
                if (i == index) { // Just inserted
                    vector_copy(c, v, invq[i]);
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

// Changes y
static void ramp(size_t c, double y[c]) {
    for (size_t i = 0; i < c; ++i) {
        if (y[i] < 0.0) { // EPS?
            y[i] = 0.0;
        }
    }
}

static void compute_u(size_t m, size_t n, size_t c, double neg_invh_f[m][n], double x[n], double neg_invh_gt[m][c], double y[c], double temp[m], double u[m]) {
    matrix_vector_product(m, n, neg_invh_f, x, u);
    ramp(c, y);
    matrix_vector_product(m, c, neg_invh_gt, y, temp); // **Sparsity
    vector_sum(m, u, temp, u);
}

// This only represents one iteration of algorithm2, i.e. algorithm2 without the lines with "while"
// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
// Note that y is modified
void algorithm2(size_t c, size_t n, size_t m, double neg_g_invh_gt[c][c], double neg_s[c][n], double neg_w[c], double neg_invh_f[m][n], double neg_invh_gt[m][c], double x[n], double invq[c][c], uint8_t a_set[c], double y[c], double v[c], double temp1[c], double temp2[m], double u[m]) {
    for (size_t i = 0; i < c; ++i) { //Later: call clear set function
        a_set[i] = 0;
    }
    matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, invq, a_set, neg_g_invh_gt, v, temp1, y);
    compute_u(m, n, c, neg_invh_f, x, neg_invh_gt, y, temp2, u);
}
