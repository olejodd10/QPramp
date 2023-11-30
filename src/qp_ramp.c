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
static size_t most_negative_index(size_t c, const iterable_set_t* a_set, double *y) {
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
static size_t most_positive_index(size_t c, const iterable_set_t* a_set, double *y) {
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

static void compute_v(size_t c, const double *invq, const iterable_set_t* a_set, size_t index, double q0, const double *neg_g_invh_gt, double *v) {
    // Compute matrix vector product
    // Sparse part
    for (size_t i = 0; i < c; ++i) {
        v[i] = set_contains(a_set, i) ? 0.0 : (i == index ? neg_g_invh_gt[c*index+i] + 1.0 : neg_g_invh_gt[c*index+i]); // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (i == index) {
            add_scaled_vector(c, v, &invq[c*i], neg_g_invh_gt[c*index+i] + 1.0, v);
        } else {
            add_scaled_vector(c, v, &invq[c*i], neg_g_invh_gt[c*index+i], v);
        }
    }
    // At this point v is defined as in the paper
    double qdiv = q0+v[index];
    if (infeasiblity_warning_enabled && 
        (fabs(qdiv) < infeasibility_warning_min || 
         fabs(qdiv) > infeasibility_warning_max)) {
        printf("WARNING: Value %e exceeds infeasibility warning limits.\n", qdiv);
    }
    scale_vector(c, v, -1.0/qdiv, v);
}

// Returns c if none found
static size_t rank_2_update_removal_index(size_t c, size_t i, const double *invq, const iterable_set_t* a_set, const double *neg_g_invh_gt, const double *y) {
    double max = 0.0; // Overwritten immediately
    size_t index = c;
    for (size_t j = set_first(a_set); j != set_end(a_set); j = set_next(a_set, j)) {
        double divisor = 0.0;
        for (size_t k = set_first(a_set); k != set_end(a_set); k = set_next(a_set, k)) {
            // Note that since j is in a_set, we only have to consider "active columns" of invq
            // Also note that the order of indices for neg_g_invh_gt doesn't matter since it's symmetric
            divisor += invq[c*k+j] * neg_g_invh_gt[c*i+k];
        }
        if ((divisor < -QP_RAMP_EPS) && (y[j]/divisor > max || index == c)) {
            max = y[j]/divisor;
            index = j;
        }
    }
    return index;
}

static void update_invq(size_t c, size_t index, const iterable_set_t* a_set, const double *v, double *invq) {
    for (size_t i = set_first(a_set); i != set_end(a_set); i = set_next(a_set, i)) {
        if (i == index) { // Just inserted
            memcpy(&invq[c*i], v, c*sizeof(double));
            invq[c*i+i] += 1.0; // Pretend there was a unit vector in the column to start with
        } else {
            add_scaled_vector(c, &invq[c*i], v, invq[c*i+index], &invq[c*i]);
        }
    }
}

static void update_y(size_t c, size_t index, const double *v, double *y) {
    add_scaled_vector(c, y, v, y[index], y);
}

// Terminal results are written to and available from a_set, invq and y
// v can be anything - call it temp?
static void algorithm1(size_t c, size_t p, double *invq, iterable_set_t* a_set, const double *neg_g_invh_gt, double *v, double *y) {
    while (1) {
        size_t index = most_negative_index(c, a_set, y);
        if (index != c) {
            compute_v(c, invq, a_set, index, 1.0, neg_g_invh_gt, v);
            set_remove(a_set, index);
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
                    printf("WARNING: Unable to perform rank two update.\n");
                } else {
                    compute_v(c, invq, a_set, index2, 1.0, neg_g_invh_gt, v);
                    set_remove(a_set, index2);
                    update_y(c, index2, v, y);
                    update_invq(c, index2, a_set, v, invq);
                }
            }

            compute_v(c, invq, a_set, index, -1.0, neg_g_invh_gt, v);
            set_insert(a_set, index);
            update_y(c, index, v, y);
            update_invq(c, index, a_set, v, invq);
        }
    }
}

static void compute_u(size_t m, size_t n, size_t c, const double *neg_invh_f, const double *x, const double *neg_g_invh, const double *y, double *u) {
    matrix_vector_product(m, n, neg_invh_f, x, u);
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > QP_RAMP_EPS) {
            add_scaled_vector(m, u, &neg_g_invh[m*i], y[i], u);
        }
    }
}

static void compute_z(size_t c, size_t p, const double *neg_g_invh, const double *y, double *z) {
    memset(z, 0, sizeof(double)*p);
    for (size_t i = 0; i < c; ++i) {
        if (y[i] > QP_RAMP_EPS) {
            add_scaled_vector(p, z, &neg_g_invh[p*i], y[i], z);
        }
    }
}

static void solve_y(size_t c, size_t n, size_t p, const double *neg_g_invh_gt, const double *neg_s, const double *neg_w, const double *x, double *y) {
    set_clear(&a_set);
    matrix_vector_product(c, n, neg_s, x, y);
    vector_sum(c, y, neg_w, y);
    algorithm1(c, p, invq, &a_set, neg_g_invh_gt, v, y);
}

void qp_ramp_solve(size_t c, size_t n, size_t p, const double *neg_g_invh_gt, const double *neg_s, const double *neg_w, const double *neg_g_invh, const double *x, double *z) {
    solve_y(c, n, p, neg_g_invh_gt, neg_s, neg_w, x, y);
    compute_z(c, p, neg_g_invh, y, z);
}

// This only represents one iteration of algorithm2, i.e. algorithm2 without the lines with "while"
// Typically you will have an MPC loop where a measured/estimated x is given as an input and then u is computed and applied
// Note that y is modified
void qp_ramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double *neg_g_invh_gt, const double *neg_s, const double *neg_w, const double *neg_invh_f, const double *neg_g_invh, const double *x, double *u) {
    solve_y(c, n, p, neg_g_invh_gt, neg_s, neg_w, x, y);
    compute_u(m, n, c, neg_invh_f, x, neg_g_invh, y, u);
}
