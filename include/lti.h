#ifndef LTI_H
#define LTI_H

#include <stddef.h>

/**
 * @brief Calculate the next state of an LTI system.
 *
 * @param[in] n Number of system states
 * @param[in] p Number of system inputs
 * @param[in] a Flattened n by n array representing system matrix
 * @param[in] x Length n array representing current system state
 * @param[in] b Flattened n by p array representing input matrix
 * @param[in] u Length p array representing system inputs
 * @param[out] res Length n array for next system state
 * @warning x and res can not be the same arrays
 */
void lti_simulate(size_t n, size_t p, const double a[n][n], const double x[n], const double b[n][p], const double u[p], double res[n]);

#endif
