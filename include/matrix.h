#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "vector.h"

// Do not use any of the functions to do operations on parts of a matrix, as the layout in memory may be unexpected and give errors.
// It seems that the array double mat[m][n] is laid out as a one-dimensional array in memory. [i][j] expands to [m*i+n] it seems.
// For example; If calling print_matrix(3,3,mat) for a matrix that is really 10x10, only values from the first row will be printed. 
// A printing function that only prints a window of a larger matrix may be implemented.

void print_matrix(size_t m, size_t n, const double *mat);

void outer_product(size_t n, const double *v1, const double *v2, double *res);

void matrix_vector_product(size_t m, size_t n, const double *mat, const double *vec, double *res);

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, const double *m1, const double *m2t, double *res);

// vec == res er lov
void negate_matrix(size_t m, size_t n, const double *mat, double *res);

// mat == res er lov
void scale_matrix(size_t m, size_t n, const double *mat, double c, double *res);

// res == m1 og res == m2 er lov
void matrix_sum(size_t m, size_t n, const double *m1, const double *m2, double *res);

// res == mat åpenbart ikke lov
// Slow!
void transpose(size_t m, size_t n, const double *mat, double *res);

uint8_t matrix_eq(size_t m, size_t n, const double *m1, const double *m2, double eps);

#endif
