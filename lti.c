#include "lti.h"

void simulate_affine(size_t n, double a[n][n], double x[n], double res[n]) {
    matrix_vector_product(n, n, a, x, res);
}

void simulate(size_t n, size_t p, double a[n][n], double x[n], double b[n][p], double u[p], double temp[n], double res[n]) {
    simulate_affine(n, a, x, temp);
    matrix_vector_product(n, p, b, u, res);
    vector_sum(n, temp, res, res);
}

void measure_affine(size_t q, size_t n, double c[q][n], double x[n], double res[q]) {
    matrix_vector_product(q, n, c, x, res);
}

void measure(size_t q, size_t n, size_t p, double c[q][n], double x[n], double d[q][p], double u[p], double temp[q], double res[q]) {
    measure_affine(q, n, c, x, temp);
    matrix_vector_product(q, p, d, u, res);
    vector_sum(q, temp, res, res);
}
