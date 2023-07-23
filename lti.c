#include "lti.h"

void simulate_affine(size_t n, const double a[n][n], const double x[n], double res[n]) {
    matrix_vector_product(n, n, a, x, res);
}

void simulate(size_t n, size_t p, const double a[n][n], const double x[n], const double b[n][p], const double u[p], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = inner_product(n, a[i], x) + inner_product(p, b[i], u);
    }
}

void measure_affine(size_t q, size_t n, const double c[q][n], const double x[n], double res[q]) {
    matrix_vector_product(q, n, c, x, res);
}

void measure(size_t q, size_t n, size_t p, const double c[q][n], const double x[n], const double d[q][p], const double u[p], double res[q]) {
    for (size_t i = 0; i < q; ++i) {
        res[i] = inner_product(n, c[i], x) + inner_product(p, d[i], u);
    }
}
