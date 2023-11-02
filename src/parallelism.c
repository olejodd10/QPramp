#include "parallelism.h"

struct inner_product_args {
    size_t n;
    const double *mat;
    const double *vec;
    double *res;
    size_t m_start;
    size_t m_end;
};

static struct inner_product_args args[NUM_THREADS];

static void *inner_product_func(void *args) {
    struct inner_product_args *args_p = (struct inner_product_args*)args;
    for (size_t i = args_p->m_start; i < args_p->m_end; ++i) {
        double result = 0.0;
        for (size_t j = 0; j < args_p->n; ++j) {
            result += args_p->mat[i*args_p->n+j] * args_p->vec[j];
        }
        args_p->res[i] = result;
    }
    return NULL;
}

void parallel_matrix_vector_product(size_t m, size_t n, const double mat[m][n], const double vec[n], double res[m]) {
    size_t chunk_size = m / NUM_THREADS;
    for (size_t i = 0 ; i < NUM_THREADS; ++i) {
    	args[i].n = n;
    	args[i].mat = (const double*)mat;
    	args[i].vec = vec;
    	args[i].res = res;
        args[i].m_start = i*chunk_size;
        args[i].m_end = i == NUM_THREADS-1 ? m : (i+1)*chunk_size;
        thread_pool_assign(i, inner_product_func, &args[i]);
    }
    thread_pool_execute();
}
