#include "iterable_set.h"

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

void set_init(iterable_set_t* set, size_t capacity) {
    set->capacity = capacity;
    set->elements = (uint8_t*)malloc(capacity*sizeof(uint8_t));
    set->next = (size_t*)malloc(capacity*sizeof(size_t));
    set->prev = (size_t*)malloc(capacity*sizeof(size_t));
    memset(set->elements, 0, capacity);
    set->first = capacity;
    set->last = capacity;
    set->size = 0;
}

void set_destroy(iterable_set_t* set) {
	free(set->elements);
	free(set->next);
	free(set->prev);
}

void set_clear(iterable_set_t* set) {
    for (size_t i = set_first(set); i != set_end(set); i = set_next(set,i)) {
        set->elements[i] = 0; // Simplified removal
    }
    set->first = set->capacity;
    set->last = set->capacity;
    set->size = 0;
}

void set_insert(iterable_set_t* set, size_t element) {
    set->elements[element] = 1;
    if (set->last == set->capacity) { // First element
        set->first = element;
        set->prev[element] = set->capacity; // Can drop this unless iterating backwards
    } else {
        set->next[set->last] = element;
        set->prev[element] = set->last;
    }
    set->next[element] = set->capacity;
    set->last = element;
    set->size++;
}

void set_remove(iterable_set_t* set, size_t element) {
    // Don't need to update prev and next for element because they only need to be correct for elements currently in the set
    set->elements[element] = 0;
    if (set->first == element) {
        set->first = set->next[element];
        if (set->last == element) {
            set->last = set->capacity;
        } else {
            set->prev[set->next[element]] = set->capacity; // Can drop this unless iterating backwards
        }
    } else if (set->last == element) {
        set->last = set->prev[element];
        set->next[set->prev[element]] = set->capacity;
    } else {
        set->next[set->prev[element]] = set->next[element];
        set->prev[set->next[element]] = set->prev[element];
    }
    set->size--;
}

size_t set_size(const iterable_set_t* set) {
    return set->size;
}

uint8_t set_contains(const iterable_set_t* set, size_t element) {
    return set->elements[element];
}

size_t set_first(const iterable_set_t* set) {
    return set->first;
}

// Iterates through active elements in order of insertion
size_t set_next(const iterable_set_t* set, size_t element) {
    return set->next[element];
}

size_t set_end(const iterable_set_t* set) {
    return set->capacity;
}
