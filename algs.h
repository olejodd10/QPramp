#ifndef ALGS_H
#define ALGS_H

#include <stdint.h>
#include <stddef.h>
#include "matrix.h"
#include "vector.h"

void algorithm2(size_t c, size_t n, size_t m, double neg_g_invh_gt[c][c], double neg_s[c][n], double neg_w[c], double neg_invh_f[m][n], double neg_g_invh[c][m], double x[n], double invq[c][c], uint8_t a_set[c], double y[c], double v[c], double u[m]);

#endif
