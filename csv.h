#ifndef CSV_H
#define CSV_H

#include <stdio.h>
#include <stddef.h>

// Returns error if file size is incompatible or if any reads fail
// Returns error if row vector is given, although the implementation would handle it
int parse_vector_csv(const char* path, size_t n, double res[n]);

// Returns error if file size is incompatible or if any reads fail
int parse_matrix_csv(const char* path, size_t m, size_t n, double res[m][n]);

// Writes as a column vector
ssize_t save_vector_csv(const char* path, size_t n, const double vec[n]);

ssize_t save_matrix_csv(const char* path, size_t m, size_t n, const double mat[m][n]);

#endif
