#include "matrix.h"

void print_matrix(size_t m, size_t n, const double mat[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        print_vector(n, mat[i]);
    }
}

void outer_product(size_t n, const double v1[n], const double v2[n], double res[n][n]) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[i][j] = v1[i]*v2[j];
        }
    }
}

void matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]) {
    for (size_t i = 0; i < m; ++i) {
        res[i] = inner_product(n, mat[i], vec);
    }
}

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, const double m1[m][n], const double m2t[p][n], double res[m][p]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < p; ++j) {
            res[i][j] = inner_product(n, m1[i], m2t[j]);
        }
    }
}

// vec == res er lov
void negate_matrix(size_t m, size_t n, const double mat[m][n], double res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        negate_vector(n, mat[i], res[i]);
    }
}

// mat == res er lov
void scale_matrix(size_t m, size_t n, const double mat[m][n], double c, double res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        scale_vector(n, mat[i], c, res[i]);
    }
}

// res == m1 og res == m2 er lov
void matrix_sum(size_t m, size_t n, const double m1[m][n], const double m2[m][n], double res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        vector_sum(n, m1[i], m2[i], res[i]);
    }
}

// res == mat åpenbart ikke lov
// Slow!
void transpose(size_t m, size_t n, const double mat[m][n], double res[n][m]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[j][i] = mat[i][j];
        }
    }
}

uint8_t matrix_eq(size_t m, size_t n, const double m1[m][n], const double m2[m][n], double eps) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(m1[i][j] - m2[i][j]) > eps) {
                return 0;
            }
        }
    }
    return 1;
}
