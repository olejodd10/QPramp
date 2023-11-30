#ifndef ALGS_H
#define ALGS_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "iterable_set.h"

// Allocates memory for solver
void qp_ramp_init(size_t c);

// Deallocates memory for solver
void qp_ramp_cleanup(void);

// Note that checking this slows down the algorithm, so consider checking output instead
void qp_ramp_enable_infeasibility_warning(double min, double max);

void qp_ramp_solve(size_t c, size_t n, size_t p, const double *neg_g_invh_gt, const double *neg_s, const double *neg_w, const double *neg_g_invh, const double *x, double *z);

void qp_ramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double *neg_g_invh_gt, const double *neg_s, const double *neg_w, const double *neg_invh_f, const double *neg_g_invh, const double *x, double *u);

#endif
