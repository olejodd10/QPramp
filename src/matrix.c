#include "matrix.h"

#include <stddef.h>

#include "vector.h"

void matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]) {
    for (size_t i = 0; i < m; ++i) {
        res[i] = vector_inner_product(n, mat[i], vec);
    }
}

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, const double m1[m][n], const double m2t[p][n], double res[m][p]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < p; ++j) {
            res[i][j] = vector_inner_product(n, m1[i], m2t[j]);
        }
    }
}

// vec == res er lov
void matrix_negate(size_t m, size_t n, const double mat[m][n], double res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        vector_negate(n, mat[i], res[i]);
    }
}

// res == mat åpenbart ikke lov
// Slow!
void matrix_transpose(size_t m, size_t n, const double mat[m][n], double res[n][m]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[j][i] = mat[i][j];
        }
    }
}
