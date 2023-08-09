#include "iterable_set.h"

void set_init(iterable_set_t* set) {
    memset(set->elements, 0, set->capacity);
    set->first = -1;
    set->last = -1;
}

void set_clear(iterable_set_t* set) {
    for (ssize_t i = set_first(set); i != -1; i = set_next(set,i)) {
        set->elements[i] = 0; // Simplified removal
    }
    set->first = -1;
    set->last = -1;
}

void set_insert(iterable_set_t* set, size_t element) {
    set->elements[element] = 1;
    if (set->last == -1) { // First element
        set->first = element;
        set->prev[element] = -1; // Can drop this unless iterating backwards
    } else {
        set->next[set->last] = element;
        set->prev[element] = set->last;
    }
    set->next[element] = -1;
    set->last = element;
}

void set_remove(iterable_set_t* set, size_t element) {
    // Don't need to update prev and next for element because they only need to be correct for elements currently in the set
    set->elements[element] = 0;
    if (set->first == element) {
        set->first = set->next[element];
        if (set->last == element) {
            set->last = -1;
        } else {
            set->prev[set->next[element]] = -1; // Can drop this unless iterating backwards
        }
    } else if (set->last == element) {
        set->last = set->prev[element];
        set->next[set->prev[element]] = -1;
    } else {
        set->next[set->prev[element]] = set->next[element];
        set->prev[set->next[element]] = set->prev[element];
    }
}

uint8_t set_contains(const iterable_set_t* set, size_t element) {
    return set->elements[element];
}

ssize_t set_first(const iterable_set_t* set) {
    return set->first;
}

// Iterates through active elements in order of insertion
ssize_t set_next(const iterable_set_t* set, size_t element) {
    return set->next[element];
}
