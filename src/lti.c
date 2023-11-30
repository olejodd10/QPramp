#include "lti.h"

void simulate_affine(size_t n, const double *a, const double *x, double *res) {
    matrix_vector_product(n, n, a, x, res);
}

void simulate(size_t n, size_t p, const double *a, const double *x, const double *b, const double *u, double *res) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = inner_product(n, &a[n*i], x) + inner_product(p, &b[p*i], u);
    }
}

void measure_affine(size_t q, size_t n, const double *c, const double *x, double *res) {
    matrix_vector_product(q, n, c, x, res);
}

void measure(size_t q, size_t n, size_t p, const double *c, const double *x, const double *d, const double *u, double *res) {
    for (size_t i = 0; i < q; ++i) {
        res[i] = inner_product(n, &c[n*i], x) + inner_product(p, &d[p*i], u);
    }
}
