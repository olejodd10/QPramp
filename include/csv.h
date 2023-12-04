#ifndef CSV_H
#define CSV_H

#include <stdio.h>
#include <stddef.h>

ssize_t csv_parse_matrix_height(const char* path); 

ssize_t csv_parse_matrix_width(const char* path); 

// Returns error if file size is incompatible or if any reads fail
// Returns error if row vector is given, although the implementation would handle it
int csv_parse_vector(const char* path, size_t n, double res[n]); 

// Returns error if file size is incompatible or if any reads fail
int csv_parse_matrix(const char* path, size_t m, size_t n, double res[m][n]); 

// Writes as a column vector
ssize_t csv_save_vector(const char* path, size_t n, const double vec[n]); 

ssize_t csv_save_matrix(const char* path, size_t m, size_t n, const double mat[m][n]); 

#endif
