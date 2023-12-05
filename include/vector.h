#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>

double vector_inner_product(size_t n, const double v1[n], const double v2[n]);

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void vector_negate(size_t n, const double vec[n], double res[n]);

// vec == res er lov
void vector_scale(size_t n, const double vec[n], double c, double res[n]);

// res == v1 og res == v2 er lov
void vector_add(size_t n, const double v1[n], const double v2[n], double res[n]);

// res == v1 og res == v2 er lov
void vector_add_scaled(size_t n, const double v1[n], const double v2[n], double c, double res[n]);

#endif
