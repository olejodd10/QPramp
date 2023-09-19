#include "thread_pool.h"

typedef struct {
    void *(*func)(void*);
    void *arg;
} task_data_t;

#define NUM_WORKERS (NUM_THREADS-1)

static size_t worker_ids[NUM_WORKERS];
static pthread_t worker_threads[NUM_WORKERS];
static pthread_barrier_t barrier;
static task_data_t task_data[NUM_THREADS];

static void* worker_thread(void* arg) {
    size_t wid = *(size_t*)arg;
    while(1) {
        pthread_barrier_wait(&barrier);
        task_data_t *assigned_task = &task_data[wid];
        assigned_task->func(assigned_task->arg);
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

void thread_pool_init() {
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    for (size_t i = 0; i < NUM_WORKERS; ++i) {
        worker_ids[i] = i;
        pthread_create(&worker_threads[i], NULL, worker_thread, &worker_ids[i]);
    }
}

void thread_pool_assign(size_t wid, void *(*func)(void*), void *arg) {
    task_data[wid] = (task_data_t) { .func = func, .arg = arg };
}

void thread_pool_execute() {
    pthread_barrier_wait(&barrier);
    task_data_t *assigned_task = &task_data[NUM_THREADS-1]; // Master thread takes last task
    assigned_task->func(assigned_task->arg);
    pthread_barrier_wait(&barrier);
}
