#ifndef ALGS_H
#define ALGS_H

#include <stdint.h>
#include "matrix.h"
#include "vector.h"

void algorithm2(uint16_t c, uint16_t n, uint16_t p, double neg_g_invh_gt[c][c], double neg_s[c][n], double neg_w[c], double neg_invh_f[p][n], double neg_invh_gt[p][c], double x[n], double invq[c][c], uint8_t a_set[c], double y[c], double v[c], double temp1[c], double temp2[c][c], double temp3[p], double u[p]);

#endif
