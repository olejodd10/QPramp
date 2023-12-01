#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <pthread.h>
#include <stddef.h>

#define NUM_THREADS 8

// Call once before scheduling tasks
void thread_pool_init();

// Schedule tasks to be run
void thread_pool_assign(size_t tid, void *(*func)(void*), void *arg);

// Run scheduled tasks and block until finished
void thread_pool_execute();

#endif
