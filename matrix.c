#include "matrix.h"

void print_matrix(uint16_t m, uint16_t n, double mat[m][n]) {
    for (int i = 0; i < m; ++i) {
        print_vector(n, mat[i]);
    }
}

void outer_product(uint16_t m, uint16_t n, double v1[m], double v2[n], double res[m][n]) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i][j] = v1[i]*v2[j];
        }
    }
}

void matrix_vector_product(uint16_t m, uint16_t n, double mat[m][n], double vec[n], double res[m]) {
    for (int i = 0; i < m; ++i) {
        res[i] = inner_product(n, mat[i], vec);
    }
}

// Merk dimensjonene til m2!! Column major
void matrix_product(uint16_t m, uint16_t n, uint16_t p, double m1[m][n], double m2[p][n], double res[m][p]) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            res[i][j] = inner_product(n, m1[i], m2[j]);
        }
    }
}

// vec == res er lov
void negate_matrix(uint16_t m, uint16_t n, double mat[m][n], double res[m][n]) {
    for (int i = 0; i < m; ++i) {
        negate_vector(n, mat[i], res[i]);
    }
}

// mat == res er lov
void scale_matrix(uint16_t m, uint16_t n, double mat[m][n], double c, double res[m][n]) {
    for (int i = 0; i < m; ++i) {
        scale_vector(n, mat[i], c, res[i]);
    }
}

// res == m1 og res == m2 er lov
void matrix_sum(uint16_t m, uint16_t n, double m1[m][n], double m2[m][n], double res[m][n]) {
    for (int i = 0; i < m; ++i) {
        vector_sum(n, m1[i], m2[i], res[i]);
    }
}

// res == mat er lov men lol
void matrix_copy(uint16_t m, uint16_t n, double mat[m][n], double res[m][n]) {
    for (int i = 0; i < m; ++i) {
        vector_copy(n, mat[i], res[i]);
    }
}

void eye(uint16_t n, double res[n][n]) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                res[i][j] = 1.0;
            } else {
                res[i][j] = 0.0;
            }
        }
    }
}

// res == mat åpenbart ikke lov
// Slow!
void transpose(uint16_t m, uint16_t n, double mat[m][n], double res[n][m]) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            res[j][i] = mat[i][j];
        }
    }
}
