#include "osqp_init.h"

size_t osqp_init_num_nonzero(size_t m, size_t n, const double mat[m][n]) {
    size_t result = 0;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(mat[i][j]) > EPS) {
                ++result;
            }
        }
    }
    return result;
}

size_t osqp_init_num_nonzero_symmetric(size_t n, const double mat[n][n]) {
    size_t result = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            if (fabs(mat[i][j]) > EPS) {
                ++result;
            }
        }
    }
    return result;
}

// mat is row major, but the result is CSC and not CSR
void osqp_init_populate_csc(size_t m, size_t n, size_t nnz, const double mat[m][n], c_float values[nnz], c_int indices[nnz], c_int ptrs[n+1]) {
    size_t num_elements = 0;
    for (size_t j = 0; j < n; ++j) {
        ptrs[j] = num_elements;
        for (size_t i = 0; i < m; ++i) {
            if (fabs(mat[i][j]) > EPS) {
                values[num_elements] = mat[i][j];
                indices[num_elements] = i; 
                ++num_elements;
            }
        }
    }
    if (nnz != num_elements) {
        printf("Incorrect number of nonzero elements was entered: %ld =/= %ld\n", nnz, num_elements);
    }
    ptrs[n] = nnz;
}

void osqp_init_populate_upper_triangular_csc(size_t n, size_t nnz, const double mat[n][n], c_float values[nnz], c_int indices[nnz], c_int ptrs[n+1]) {
    size_t num_elements = 0;
    for (size_t j = 0; j < n; ++j) {
        ptrs[j] = num_elements;
        for (size_t i = 0; i <= j; ++i) {
            if (fabs(mat[i][j]) > EPS) {
                values[num_elements] = mat[i][j];
                indices[num_elements] = i; 
                ++num_elements;
            }
        }
    }
    if (nnz != num_elements) {
        printf("Incorrect number of nonzero elements was entered: %ld =/= %ld\n", nnz, num_elements);
    }
    ptrs[n] = nnz;
}

