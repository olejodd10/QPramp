#include "qpramp.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "vector.h"
#include "iterable_set.h"
#include "indexed_vectors.h"

#define QPRAMP_EPS 1e-8

static uint8_t infeasiblity_warning_enabled = 0;
static double infeasibility_warning_min;
static double infeasibility_warning_max;

static double *y = NULL;
static double *v = NULL;
static iterable_set_t a_set;
static indexed_vectors_t invq;

void qpramp_init(size_t c, size_t p) {
	y = (double*)malloc(c*sizeof(double));
	v = (double*)malloc(c*sizeof(double));
    set_init(&a_set, c);
    indexed_vectors_init(&invq, p, c, c);
}

void qpramp_cleanup(void) {
	free(y);
	free(v);
    set_destroy(&a_set);
    indexed_vectors_destroy(&invq);
}

void qpramp_enable_infeasibility_error(double min, double max) {
    infeasiblity_warning_enabled = 1;
    infeasibility_warning_min = min;
    infeasibility_warning_max = max;
}

// Returns c if none found
static inline size_t most_negative_index(size_t c, const iterable_set_t* a_set, double y[c]) {
    double min = -QPRAMP_EPS;
    size_t index = c; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (y[i] < min) {
            min = y[i];
            index = i;
        }    
    }
    return index;
}

// Returns c if none found
static inline size_t most_positive_index(size_t c, const iterable_set_t* a_set, double y[c]) {
    double max = QPRAMP_EPS;
    size_t index = c; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > max && !set_contains(a_set, i)) {
            max = y[i];
            index = i;
        }    
    }
    return index;
}

static int compute_v(size_t c, const indexed_vectors_t *invq, const iterable_set_t* a_set, size_t index, double q0, const double neg_g_invh_gt[c][c], double v[c]) {
    // Compute matrix vector product
    // Sparse part
    for (size_t i = 0; i < c; ++i) {
        v[i] = set_contains(a_set, i) ? 0.0 : (i == index ? neg_g_invh_gt[index][i] + 1.0 : neg_g_invh_gt[index][i]); // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (i == index) {
            vector_add_scaled(c, v, indexed_vectors_get(invq, i), neg_g_invh_gt[index][i] + 1.0, v);
        } else {
            vector_add_scaled(c, v, indexed_vectors_get(invq, i), neg_g_invh_gt[index][i], v);
        }
    }
    // At this point v is defined as in the paper
    double qdiv = q0+v[index];
    if (infeasiblity_warning_enabled && 
        (fabs(qdiv) < infeasibility_warning_min || 
         fabs(qdiv) > infeasibility_warning_max)) {
        return QPRAMP_ERROR_INFEASIBLE;
    }
    vector_scale(c, v, -1.0/qdiv, v);
    return 0;
}

// Returns c if none found
static inline size_t rank_2_update_removal_index(size_t c, size_t i, const indexed_vectors_t *invq, const iterable_set_t* a_set, const double neg_g_invh_gt[c][c], const double y[c]) {
    double max = 0.0;
    size_t index = c;
    for (size_t j = set_first(a_set); j != set_end(a_set); j = set_next(a_set, j)) {
        double divisor = 0.0;
        for (size_t k = set_first(a_set); k != set_end(a_set); k = set_next(a_set, k)) {
            // Note that the order of indices for neg_g_invh_gt doesn't matter since it's symmetric
            divisor += indexed_vectors_get(invq, k)[j] * neg_g_invh_gt[i][k];
        }
        if ((divisor < -QPRAMP_EPS) && (y[j]/divisor > max || index == c)) {
            max = y[j]/divisor;
            index = j;
        }
    }
    return index;
}

static inline void update_invq(size_t c, size_t index, const iterable_set_t* a_set, const double v[c], indexed_vectors_t *invq) {
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        double* invqi = indexed_vectors_get_mut(invq, i);
        vector_add_scaled(c, invqi, v, invqi[index], invqi);
    }
}

static inline void update_y(size_t c, size_t index, const double v[c], double y[c]) {
    vector_add_scaled(c, y, v, y[index], y);
}

static int algorithm1(size_t c, size_t p, indexed_vectors_t *invq, iterable_set_t* a_set, const double neg_g_invh_gt[c][c], double v[c], double y[c]) {
    while (1) {
        size_t index = most_negative_index(c, a_set, y);
        if (index != c) {
            int err = compute_v(c, invq, a_set, index, 1.0, neg_g_invh_gt, v);
            if (err) {
                return err;
            }
            set_remove(a_set, index);
            indexed_vectors_remove(invq, index);
            update_y(c, index, v, y);
            update_invq(c, index, a_set, v, invq);
        } else {
            index = most_positive_index(c, a_set, y);
            if (index == c) {
                break;
            }

            if (set_size(a_set) == p) {
                size_t index2 = rank_2_update_removal_index(c, index, invq, a_set, neg_g_invh_gt, y);
                if (index2 == c) {
                    return QPRAMP_ERROR_RANK_2_UPDATE;
                } else {
                    int err = compute_v(c, invq, a_set, index2, 1.0, neg_g_invh_gt, v);
                    if (err) {
                        return err;
                    }
                    set_remove(a_set, index2);
                    indexed_vectors_remove(invq, index2);
                    update_y(c, index2, v, y);
                    update_invq(c, index2, a_set, v, invq);
                }
            }

            int err = compute_v(c, invq, a_set, index, -1.0, neg_g_invh_gt, v);
            if (err) {
                return err;
            }
            update_y(c, index, v, y);
            update_invq(c, index, a_set, v, invq);
            set_insert(a_set, index);
            indexed_vectors_insert(invq, index, v);
            indexed_vectors_get_mut(invq, index)[index] += 1.0; // Pretend there was a unit vector in the column to start with
        }
    }
    return 0;
}

static void compute_u(size_t m, size_t n, size_t c, const iterable_set_t *a_set, const double neg_invh_f[m][n], const double x[n], const double neg_g_invh[c][m], const double y[c], double u[m]) {
    matrix_vector_product(m, n, neg_invh_f, x, u);
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        vector_add_scaled(m, u, neg_g_invh[i], y[i], u);
    }
}

static void compute_z(size_t c, size_t p, const iterable_set_t *a_set, const double neg_g_invh[c][p], const double y[c], double z[p]) {
    memset(z, 0, sizeof(double)*p);
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        vector_add_scaled(p, z, neg_g_invh[i], y[i], z);
    }
}

static inline void initialize_y(size_t c, size_t n, const indexed_vectors_t *invq, const iterable_set_t *a_set, const double neg_s[c][n], const double neg_w[c], const double x[n], double v[c], double y[c]) {
    matrix_vector_product(c, n, neg_s, x, v); // Use v temporarily; it will be overwritten later
    vector_add(c, v, neg_w, v);
    memcpy(y, v, c*sizeof(double));
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        vector_add_scaled(c, y, indexed_vectors_get(invq, i), v[i], y);
        y[i] -= v[i];
    }
}

static int solve_y(size_t c, size_t n, size_t p, indexed_vectors_t *invq, iterable_set_t *a_set, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double x[n], double v[c], double y[c]) {
    initialize_y(c, n, invq, a_set, neg_s, neg_w, x, v, y);
    return algorithm1(c, p, invq, a_set, neg_g_invh_gt, v, y);
}

int qpramp_solve(size_t c, size_t n, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_g_invh[c][p], const double x[n], double z[p]) {
    int err = solve_y(c, n, p, &invq, &a_set, neg_g_invh_gt, neg_s, neg_w, x, v, y);
    if (err) {
        return err;
    }
    compute_z(c, p, &a_set, neg_g_invh, y, z);
    return 0;
}

// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
int qpramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_invh_f[m][n], const double neg_g_invh[c][m], const double x[n], double u[m]) {
    int err = solve_y(c, n, p, &invq, &a_set, neg_g_invh_gt, neg_s, neg_w, x, v, y);
    if (err) {
        return err;
    }
    compute_u(m, n, c, &a_set, neg_invh_f, x, neg_g_invh, y, u);
    return 0;
}
