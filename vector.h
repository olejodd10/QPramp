#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stddef.h>

void print_vector(size_t n, double vec[n]); 

double inner_product(size_t n, double v1[n], double v2[n]); 

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void negate_vector(size_t n, double vec[n], double res[n]); 

// vec == res er lov
void scale_vector(size_t n, double vec[n], double c, double res[n]); 

// res == v1 og res == v2 er lov
void vector_sum(size_t n, double v1[n], double v2[n], double res[n]);

// res == vec er lov men lol
void vector_copy(size_t n, double vec[n], double res[n]); 

#endif
