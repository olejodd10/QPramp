#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stddef.h>

void print_vector(size_t n, const double *vec);

double inner_product(size_t n, const double *v1, const double *v2);

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void negate_vector(size_t n, const double *vec, double *res);

// vec == res er lov
void scale_vector(size_t n, const double *vec, double c, double *res);

// res == v1 og res == v2 er lov
void vector_sum(size_t n, const double *v1, const double *v2, double *res);

// res == v1 og res == v2 er lov
void add_scaled_vector(size_t n, const double *v1, const double *v2, double c, double *res);

#endif
