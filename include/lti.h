#ifndef LTI_H
#define LTI_H

#include <stddef.h>

// res == x IKKE lov
void lti_simulate(size_t n, size_t p, const double a[n][n], const double x[n], const double b[n][p], const double u[p], double res[n]);

#endif
