#include "lti.h"

void simulate_affine(size_t n, double a[n][n], double x[n], double res[n]) {
    matrix_vector_product(n, n, a, x, res);
}

void simulate(size_t n, size_t p, double a[n][n], double x[n], double b[n][p], double u[p], double temp[n], double res[n]) {
    simulate_affine(n, a, x, res);
    matrix_vector_product(n, p, b, u, temp);
    vector_sum(n, res, temp, res);
}

void measure_affine(size_t q, size_t n, double c[q][n], double x[n], double res[q]) {
    matrix_vector_product(q, n, c, x, res);
}

void measure(size_t q, size_t n, size_t p, double c[q][n], double x[n], double d[q][p], double u[p], double temp[q], double res[q]) {
    measure_affine(q, n, c, x, res);
    matrix_vector_product(q, p, d, u, temp);
    vector_sum(q, res, temp, res);
}
