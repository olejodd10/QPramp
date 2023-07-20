#include "timing.h"

clock_t start;

void tick() {
    start = clock();
}

int tock() {
    return (clock() - start) * 1000 / CLOCKS_PER_SEC;
}
