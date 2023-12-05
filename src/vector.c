#include "vector.h"

#include <stddef.h>

double vector_inner_product(size_t n, const double v1[n], const double v2[n]) {
    double result = 0.0;
    for (size_t i = 0; i < n; ++i) {
        result += v1[i]*v2[i];
    }
    return result;
}

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void vector_negate(size_t n, const double vec[n], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = -vec[i];
    }
}

// vec == res er lov
void vector_scale(size_t n, const double vec[n], double c, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = c*vec[i];
    }
}

// res == v1 og res == v2 er lov
void vector_add(size_t n, const double v1[n], const double v2[n], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = v1[i] + v2[i];
    }
}

// res == v1 og res == v2 er lov
void vector_add_scaled(size_t n, const double v1[n], const double v2[n], double c, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = v1[i] + c*v2[i];
    }
}
