#include "lti.h"

#include <stddef.h>

#include "vector.h"

void lti_simulate(size_t n, size_t p, const double a[n][n], const double x[n], const double b[n][p], const double u[p], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = inner_product(n, a[i], x) + inner_product(p, b[i], u);
    }
}
