#include "timing.h"
#include <time.h>

#define TIMING_CLOCK CLOCK_REALTIME

static struct timespec then, now;

void timing_print_precision() {
    clock_getres(TIMING_CLOCK, &then);
    printf("Timing precision is %ld s %ld ns\n", then.tv_sec, then.tv_nsec);
}

void timing_reset() {
    clock_gettime(TIMING_CLOCK, &then);
}

long timing_elapsed() {
    clock_gettime(TIMING_CLOCK, &now);
    return (now.tv_sec-then.tv_sec)*1000000000 + now.tv_nsec-then.tv_nsec;
}
