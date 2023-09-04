#ifndef OSQP_INIT_H
#define OSQP_INIT_H

#include <math.h>
#include <stddef.h>
#include "osqp.h"

size_t osqp_init_num_nonzero(size_t m, size_t n, const double mat[m][n]);

size_t osqp_init_num_nonzero_symmetric(size_t n, const double mat[n][n]);

// mat is row major, but the result is CSC and not CSR
void osqp_init_populate_csc(size_t m, size_t n, size_t nnz, const double mat[m][n], c_float values[nnz], c_int indices[nnz], c_int ptrs[n+1]);

void osqp_init_populate_upper_triangular_csc(size_t n, size_t nnz, const double mat[n][n], c_float values[nnz], c_int indices[nnz], c_int ptrs[n+1]);

#endif 
