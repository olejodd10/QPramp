#include "qp_ramp.h"

#define QP_RAMP_EPS 1e-8

static uint8_t infeasiblity_warning_enabled = 0;
static double infeasibility_warning_min;
static double infeasibility_warning_max;

static double *y = NULL;
static double *v = NULL;
static double *invq = NULL;
static iterable_set_t a_set;

void qp_ramp_init(size_t c) {
	y = (double*)malloc(c*sizeof(double));
	v = (double*)malloc(c*sizeof(double));
    invq = (double*)malloc(c*c*sizeof(double));
    set_init(&a_set, c);
    thread_pool_init();
}

void qp_ramp_cleanup(void) {
	free(y);
	free(v);
	free(invq);
    set_destroy(&a_set);
}

void qp_ramp_enable_infeasibility_warning(double min, double max) {
    infeasiblity_warning_enabled = 1;
    infeasibility_warning_min = min;
    infeasibility_warning_max = max;
}

// Returns c if none found
static size_t most_negative_index(size_t c, const iterable_set_t* a_set, double y[c]) {
    double min = -QP_RAMP_EPS;
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
static size_t most_positive_index(size_t c, const iterable_set_t* a_set, double y[c]) {
    double max = QP_RAMP_EPS;
    size_t index = c; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > max && !set_contains(a_set, i)) {
            max = y[i];
            index = i;
        }    
    }
    return index;
}

static void compute_v(size_t c, const double invq[c][c], const iterable_set_t* a_set, size_t index, const double neg_g_invh_gt[c][c], double v[c]) {
    // Compute matrix vector product
    // Sparse part
    for (size_t i = 0; i < c; ++i) {
        v[i] = set_contains(a_set, i) ? 0.0 : (i == index ? neg_g_invh_gt[index][i] + 1.0 : neg_g_invh_gt[index][i]); // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (i == index) {
            add_scaled_vector(c, v, invq[i], neg_g_invh_gt[index][i] + 1.0, v);
        } else {
            add_scaled_vector(c, v, invq[i], neg_g_invh_gt[index][i], v);
        }
    }
}

// Returns c if none found
static size_t rank_2_update_removal_index(size_t c, size_t i, const double invq[c][c], const iterable_set_t* a_set, const double neg_g_invh_gt[c][c], const double y[c]) {
    double max = 0.0; // Overwritten immediately
    size_t index = c;
    for (size_t j = set_first(a_set); j != set_end(a_set); j = set_next(a_set, j)) {
        double divisor = 0.0;
        for (size_t k = set_first(a_set); k != set_end(a_set); k = set_next(a_set, k)) {
            // Note that since j is in a_set, we only have to consider "active columns" of invq
            // Also note that the order of indices for neg_g_invh_gt doesn't matter since it's symmetric
            divisor += invq[k][j] * neg_g_invh_gt[i][k];
        }
        if ((divisor < -QP_RAMP_EPS) && (y[j]/divisor > max || index == c)) {
            max = y[j]/divisor;
            index = j;
        }
    }
    return index;
}

// Note that v is also changed
static void update_y_and_invq(size_t c, size_t index, double q0, const iterable_set_t* a_set, double v[c], double y[c], double invq[c][c]) {
    double qdiv = q0+v[index];
    if (infeasiblity_warning_enabled && 
        (fabs(qdiv) < infeasibility_warning_min || 
         fabs(qdiv) > infeasibility_warning_max)) {
        printf("WARNING: Value %e exceeds infeasibility warning limits.\n", qdiv);
    }
    // Update invq
    scale_vector(c, v, -1.0/qdiv, v);
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (i == index) { // Just inserted
            memcpy(invq[i], v, c*sizeof(double));
            invq[i][i] += 1.0; // Pretend there was a unit vector in the column to start with
        } else {
            add_scaled_vector(c, invq[i], v, invq[i][index], invq[i]);
        }
    }
    // Update y
    add_scaled_vector(c, y, v, y[index], y);
}

// Terminal results are written to and available from a_set, invq and y
// v can be anything - call it temp?
static void algorithm1(size_t c, size_t p, double invq[c][c], iterable_set_t* a_set, const double neg_g_invh_gt[c][c], double v[c], double y[c]) {
    while (1) {
        size_t index = most_negative_index(c, a_set, y);
        if (index != c) {
            compute_v(c, invq, a_set, index, neg_g_invh_gt, v);
            set_remove(a_set, index);
            update_y_and_invq(c, index, 1.0, a_set, v, y, invq);
        } else {
            index = most_positive_index(c, a_set, y);
            if (index == c) {
                break;
            }

            if (set_size(a_set) == p) {
                size_t index2 = rank_2_update_removal_index(c, index, invq, a_set, neg_g_invh_gt, y);
                if (index2 == c) {
                    printf("WARNING: Unable to perform rank two update.\n");
                } else {
                    compute_v(c, invq, a_set, index2, neg_g_invh_gt, v);
                    set_remove(a_set, index2);
                    update_y_and_invq(c, index2, 1.0, a_set, v, y, invq);
                }
            }

            compute_v(c, invq, a_set, index, neg_g_invh_gt, v);
            set_insert(a_set, index);
            update_y_and_invq(c, index, -1.0, a_set, v, y, invq);
        }
    }
}

static void compute_u(size_t m, size_t n, size_t c, const double neg_invh_f[m][n], const double x[n], const double neg_g_invh[c][m], const double y[c], double u[m]) {
    matrix_vector_product(m, n, neg_invh_f, x, u); // Typically few inputs, so parallelization overhead is not worth it
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > QP_RAMP_EPS) {
            add_scaled_vector(m, u, neg_g_invh[i], y[i], u);
        }
    }
}

static void compute_z(size_t c, size_t p, const double neg_g_invh[c][p], const double y[c], double z[p]) {
    memset(z, 0, sizeof(double)*p);
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > QP_RAMP_EPS) {
            add_scaled_vector(p, z, neg_g_invh[i], y[i], z);
        }
    }
}

void qp_ramp_solve(size_t c, size_t n, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_g_invh[c][p], const double x[n], double z[p]) {
    set_clear(&a_set);
    matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, p, (double(*)[])invq, &a_set, neg_g_invh_gt, v, y);
    compute_z(c, p, neg_g_invh, y, z);
}

// This only represents one iteration of algorithm2, i.e. algorithm2 without the lines with "while"
// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
// Note that y is modified
void qp_ramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_invh_f[m][n], const double neg_g_invh[c][m], const double x[n], double u[m]) {
    set_clear(&a_set);
    parallel_matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, p, (double(*)[])invq, &a_set, neg_g_invh_gt, v, y);
    compute_u(m, n, c, neg_invh_f, x, neg_g_invh, y, u);
}
