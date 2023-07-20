#ifndef PARSE_H
#define PARSE_H

#include <stdio.h>
#include <stddef.h>

// Leaves stream position unchanged
size_t parse_matrix_csv_height(FILE* f);

// Leaves stream position unchanged
size_t parse_matrix_csv_width(FILE* f);

// Leaves stream position unchanged
int assert_size_csv(FILE* f, size_t m, size_t n);

// Assumes stream is at start of file. Changes stream position.
// Returns error if file size is incompatible or if any reads fail
// Returns error if row vector is given, although the implementation would handle it
int parse_vector_csv(FILE* f, size_t n, double res[n]);

// Assumes stream is at start of file. Changes stream position.
// Returns error if file size is incompatible or if any reads fail
int parse_matrix_csv(FILE* f, size_t m, size_t n, double res[m][n]);

#endif
