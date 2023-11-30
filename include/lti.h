#ifndef LTI_H
#define LTI_H

#include <stdint.h>
#include <stddef.h>
#include "matrix.h"
#include "vector.h"

// res == x IKKE lov
void simulate_affine(size_t n, const double *a, const double *x, double *res);

// res == x IKKE lov
void simulate(size_t n, size_t p, const double *a, const double *x, const double *b, const double *u, double *res);

void measure_affine(size_t q, size_t n, const double *c, const double *x, double *res);

void measure(size_t q, size_t n, size_t p, const double *c, const double *x, const double *d, const double *u, double *res);

#endif
