#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <stdio.h>

void print_vector(uint16_t n, double vec[n]); 

double inner_product(uint16_t n, double v1[n], double v2[n]); 

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void negate_vector(uint16_t n, double vec[n], double res[n]); 

// vec == res er lov
void scale_vector(uint16_t n, double vec[n], double c, double res[n]); 

// res == v1 og res == v2 er lov
void vector_sum(uint16_t n, double v1[n], double v2[n], double res[n]);

// res == vec er lov men lol
void vector_copy(uint16_t n, double vec[n], double res[n]); 

#endif
