#ifndef PARALLELISM_H
#define PARALLELISM_H

#include <stddef.h>
#include "matrix.h"
#include "vector.h"
#include "thread_pool.h"

void parallel_matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]);

#endif
