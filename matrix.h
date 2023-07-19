#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include "vector.h"

void print_matrix(size_t m, size_t n, double mat[m][n]);

void outer_product(size_t m, size_t n, double v1[m], double v2[n], double res[m][n]);

void matrix_vector_product(size_t m, size_t n, double mat[m][n], double vec[n], double res[m]);

// Merk dimensjonene til m2!! Column major/transponerte av m2 
void matrix_product(size_t m, size_t n, size_t p, double m1[m][n], double m2t[p][n], double res[m][p]);

// vec == res er lov
void negate_matrix(size_t m, size_t n, double mat[m][n], double res[m][n]);

// mat == res er lov
void scale_matrix(size_t m, size_t n, double mat[m][n], double c, double res[m][n]);

// res == m1 og res == m2 er lov
void matrix_sum(size_t m, size_t n, double m1[m][n], double m2[m][n], double res[m][n]);

// res == mat er lov men lol
void matrix_copy(size_t m, size_t n, double mat[m][n], double res[m][n]);

void eye(size_t n, double res[n][n]);

// res == mat åpenbart ikke lov
// Slow!
void transpose(size_t m, size_t n, double mat[m][n], double res[n][m]);

#endif
