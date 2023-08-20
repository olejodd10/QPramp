#include "thread_pool.h"
#include <stdio.h>

typedef struct {
    void *(*func)(void*);
    void *arg;
} task_data_t;

static pthread_t threads[NUM_THREADS];
static pthread_barrier_t barrier;
static task_data_t task_data[NUM_THREADS];

static size_t thread_ids[NUM_THREADS];

static void* worker_thread(void* arg) {
    size_t tid = *(size_t*)arg;
    while(1) {
        pthread_barrier_wait(&barrier);
        task_data_t *assigned_task = &task_data[tid];
        assigned_task->func(assigned_task->arg);
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

void thread_pool_init() {
    pthread_barrier_init(&barrier, NULL, NUM_THREADS+1);
    for (size_t i = 0; i < NUM_THREADS; ++i) {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, worker_thread, &thread_ids[i]);
    }
}

void thread_pool_assign(size_t tid, void *(*func)(void*), void *arg) {
    task_data[tid] = (task_data_t) { .func = func, .arg = arg };
}

void thread_pool_execute() {
    pthread_barrier_wait(&barrier);
    pthread_barrier_wait(&barrier);
}
