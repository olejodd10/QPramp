#include "matrix.h"

void print_matrix(size_t m, size_t n, const double *mat) {
    for (size_t i = 0; i < m; ++i) {
        print_vector(n, &mat[n*i]);
    }
}

void outer_product(size_t n, const double *v1, const double *v2, double *res) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[n*i+j] = v1[i]*v2[j];
        }
    }
}

void matrix_vector_product(size_t m, size_t n, const double *mat, const double *vec, double *res) {
    for (size_t i = 0; i < m; ++i) {
        res[i] = inner_product(n, &mat[n*i], vec);
    }
}

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, const double *m1, const double *m2t, double *res) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < p; ++j) {
            res[p*i+j] = inner_product(n, &m1[n*i], &m2t[n*j]);
        }
    }
}

// vec == res er lov
void negate_matrix(size_t m, size_t n, const double *mat, double *res) {
    for (size_t i = 0; i < m; ++i) {
        negate_vector(n, &mat[n*i], &res[n*i]);
    }
}

// mat == res er lov
void scale_matrix(size_t m, size_t n, const double *mat, double c, double *res) {
    for (size_t i = 0; i < m; ++i) {
        scale_vector(n, &mat[n*i], c, &res[n*i]);
    }
}

// res == m1 og res == m2 er lov
void matrix_sum(size_t m, size_t n, const double *m1, const double *m2, double *res) {
    for (size_t i = 0; i < m; ++i) {
        vector_sum(n, &m1[n*i], &m2[n*i], &res[n*i]);
    }
}

// res == mat åpenbart ikke lov
// Slow!
void transpose(size_t m, size_t n, const double *mat, double *res) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[m*j+i] = mat[n*i+j];
        }
    }
}

uint8_t matrix_eq(size_t m, size_t n, const double *m1, const double *m2, double eps) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(m1[n*i+j] - m2[n*i+j]) > eps) {
                return 0;
            }
        }
    }
    return 1;
}
