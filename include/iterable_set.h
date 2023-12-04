#ifndef ITERABLE_SET_H
#define ITERABLE_SET_H

#include <stdint.h>
#include <stddef.h>

typedef struct {
    size_t capacity;
    size_t size;
    uint8_t* elements;
    size_t* next;
    size_t* prev;
    size_t first;
    size_t last;
} iterable_set_t;

void set_init(iterable_set_t* set, size_t capacity);

void set_destroy(iterable_set_t* set);

void set_clear(iterable_set_t* set);

void set_insert(iterable_set_t* set, size_t element);

void set_remove(iterable_set_t* set, size_t element);

size_t set_size(const iterable_set_t* set);

uint8_t set_contains(const iterable_set_t* set, size_t element);

size_t set_first(const iterable_set_t* set);

// Iterates through active elements in order of insertion
size_t set_next(const iterable_set_t* set, size_t element);

// Returns the value returned by set_next when called with the last element
size_t set_end(const iterable_set_t* set);

#endif
