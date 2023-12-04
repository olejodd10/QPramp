#ifndef ALGS_H
#define ALGS_H

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "iterable_set.h"
#include "indexed_vectors.h"

#define QPRAMP_ERROR_INFEASIBLE -1
#define QPRAMP_ERROR_RANK_2_UPDATE -2

// Allocates memory for solver
void qp_ramp_init(size_t c, size_t p);

// Deallocates memory for solver
void qp_ramp_cleanup(void);

// Note that checking this slows down the algorithm, so consider checking output instead
void qp_ramp_enable_infeasibility_error(double min, double max);

int qp_ramp_solve(size_t c, size_t n, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_g_invh[c][p], const double x[n], double z[p]);

int qp_ramp_solve_mpc(size_t c, size_t n, size_t m, size_t p, const double neg_g_invh_gt[c][c], const double neg_s[c][n], const double neg_w[c], const double neg_invh_f[m][n], const double neg_g_invh[c][m], const double x[n], double u[m]);

#endif
