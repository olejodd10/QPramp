#ifndef PARSE_H
#define PARSE_H

#include <stdio.h>
#include <stddef.h>

// Resets stream position, reads and then resets position to original
size_t parse_matrix_csv_height(FILE* f);

// Resets stream position, reads and then resets position to original
size_t parse_matrix_csv_width(FILE* f);

// Assumes stream is at start of vector. Does not reset.
// Note that this works with both row- and column vectors
void parse_vector_csv(FILE* f, size_t n, double res[n]);

// Assumes stream is at start of matrix. Does not reset.
void parse_matrix_csv(FILE* f, size_t m, size_t n, double res[m][n]);

#endif
