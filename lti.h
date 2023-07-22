#ifndef LTI_H
#define LTI_H

#include <stdint.h>
#include <stddef.h>
#include "matrix.h"
#include "vector.h"

// res == x IKKE lov
void simulate_affine(size_t n, double a[n][n], double x[n], double res[n]);

// res == x IKKE lov
void simulate(size_t n, size_t p, double a[n][n], double x[n], double b[n][p], double u[p], double res[n]); 

void measure_affine(size_t q, size_t n, double c[q][n], double x[n], double res[q]); 

void measure(size_t q, size_t n, size_t p, double c[q][n], double x[n], double d[q][p], double u[p], double res[q]); 

#endif
