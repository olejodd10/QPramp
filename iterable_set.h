#ifndef ITERABLE_SET_H
#define ITERABLE_SET_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    size_t capacity;
    uint8_t* elements;
    ssize_t* next;
    ssize_t* prev;
    ssize_t first;
    ssize_t last;
} iterable_set_t;

void set_init(iterable_set_t* set);

void set_clear(iterable_set_t* set);

void set_insert(iterable_set_t* set, size_t element);

void set_remove(iterable_set_t* set, size_t element);

uint8_t set_contains(const iterable_set_t* set, size_t element);

ssize_t set_first(const iterable_set_t* set);

// Iterates through active elements in order of insertion
ssize_t set_next(const iterable_set_t* set, size_t element);

#endif
