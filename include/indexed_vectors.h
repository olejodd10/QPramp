#ifndef MATRIX_COLUMNS_H
#define MATRIX_COLUMNS_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    size_t length;
    size_t capacity;
    size_t max_index;
    size_t size;
    size_t* whats_here;
    size_t* where_is;
    double* values;
} indexed_vectors_t;

void indexed_vectors_clear(indexed_vectors_t *vecs);

void indexed_vectors_init(indexed_vectors_t *vecs, size_t capacity, size_t length, size_t max_index);

void indexed_vectors_destroy(indexed_vectors_t *vecs);

static size_t first_free(const indexed_vectors_t *vecs);

double* indexed_vectors_get_mut(indexed_vectors_t *vecs, size_t index);

const double* indexed_vectors_get(const indexed_vectors_t *vecs, size_t index);

void indexed_vectors_insert(indexed_vectors_t *vecs, size_t index, const double *column);

void indexed_vectors_remove(indexed_vectors_t *vecs, size_t index);

#endif
