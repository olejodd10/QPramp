#ifndef LTI_H
#define LTI_H

#include <stdint.h>
#include "matrix.h"
#include "vector.h"

void simulate_affine(uint16_t n, double a[n][n], double x[n], double res[n]);

void simulate(uint16_t n, uint16_t p, double a[n][n], double x[n], double b[n][p], double u[p], double temp[n], double res[n]);

void measure_affine(uint16_t q, uint16_t n, double c[q][n], double x[n], double res[q]);

void measure(uint16_t q, uint16_t n, uint16_t p, double c[q][n], double x[n], double d[q][p], double u[p], double temp[q], double res[q]);

#endif
