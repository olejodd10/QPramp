#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include "vector.h"

// Do not use any of the functions to do operations on parts of a matrix, as the layout in memory may be unexpected and give errors.
// It seems that the array double mat[m][n] is laid out as a one-dimensional array in memory. [i][j] expands to [m*i+n] it seems.
// For example; If calling print_matrix(3,3,mat) for a matrix that is really 10x10, only values from the first row will be printed. 
// A printing function that only prints a window of a larger matrix may be implemented.

void print_matrix(size_t m, size_t n, const double mat[m][n]);

void outer_product(size_t n, const double v1[n], const double v2[n], double res[n][n]);

void matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]);

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, const double m1[m][n], const double m2t[p][n], double res[m][p]);

// vec == res er lov
void negate_matrix(size_t m, size_t n, const double mat[m][n], double res[m][n]);

// mat == res er lov
void scale_matrix(size_t m, size_t n, const double mat[m][n], double c, double res[m][n]);

// res == m1 og res == m2 er lov
void matrix_sum(size_t m, size_t n, const double m1[m][n], const double m2[m][n], double res[m][n]);

// res == mat åpenbart ikke lov
// Slow!
void transpose(size_t m, size_t n, const double mat[m][n], double res[n][m]);

#endif
