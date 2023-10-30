#ifndef LTI_H
#define LTI_H

#include <stdint.h>
#include <stddef.h>
#include "matrix.h"
#include "vector.h"

// res == x IKKE lov
void simulate_affine(size_t n, const double a[n][n], const double x[n], double res[n]);

// res == x IKKE lov
void simulate(size_t n, size_t p, const double a[n][n], const double x[n], const double b[n][p], const double u[p], double res[n]);

void measure_affine(size_t q, size_t n, const double c[q][n], const double x[n], double res[q]);

void measure(size_t q, size_t n, size_t p, const double c[q][n], const double x[n], const double d[q][p], const double u[p], double res[q]);

#endif
