#include "vector.h"

void print_vector(size_t n, double vec[n]) {
    for (size_t i = 0; i < n; ++i) {
        printf("%.2e ", vec[i]);
    }
    printf("\n");
}

double inner_product(size_t n, double v1[n], double v2[n]) {
    double result = 0.0;
    for (size_t i = 0; i < n; ++i) {
        result += v1[i]*v2[i];
    }
    return result;
}

// https://stackoverflow.com/questions/3412190/net-decimal-negate-vs-multiplying-by-1
// vec == res er lov
void negate_vector(size_t n, double vec[n], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = -vec[i];
    }
}

// vec == res er lov
void scale_vector(size_t n, double vec[n], double c, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = c*vec[i];
    }
}

// res == v1 og res == v2 er lov
void vector_sum(size_t n, double v1[n], double v2[n], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = v1[i] + v2[i];
    }
}

// res == vec er lov men lol
void vector_copy(size_t n, double vec[n], double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = vec[i];
    }
}

