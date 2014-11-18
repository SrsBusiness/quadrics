#include "quadric.h"
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "minicurses.h"
#include <pthread.h>
#include <assert.h>

void display_subspace(subspace *);

void single_thread_benchmark(int64_t);
void multi_thread_benchmark(int64_t, int);
void find_surface_test(int64_t);
void herp_test();

int main(int argc, char **argv) {
    multi_thread_benchmark(19, 100);
    //single_thread_benchmark(19);
    //find_surface_test(20);
    //herp_test();
}

void herp_test() {
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -19 * 19};
    vector v = {-10, 16, 2};
    is_surface(&q, &v);
}

void find_surface_test(int64_t radius) {
    subspace *s = subspace_init(-radius - 1, -radius - 1, 0, radius + 2,
            radius + 2, 1);
    quadric q = {1, 1, 0, 0, 0, 0, 0, 0, 0, -radius * radius};
    vector v = {0, 0, 0};
    vector surface;
    assert(find_surface(&q, &v, &surface)); 
    breadth_first_fill(s, &q, &surface);
    display_subspace(s);
}

void print_elem(void *elem) {
    printf("%ld\n", (int64_t)elem);
}

typedef struct _builder_args {
    subspace *s;
    quadric *q;
    uint64_t index;
    void (*func)(subspace *, const quadric *, const vector *);
} builder_args;

void *builder_thread(void *args) {
    builder_args *bargs = (builder_args *)args;
    vector v, surface;
    //v.x = _x(bargs->s, bargs->index);
    //v.y = _y(bargs->s, bargs->index);
    //v.z = _z(bargs->s, bargs->index);
    v.x = -10;
    v.y = 12;
    v.z = 0;
    int success = find_surface(bargs->q, &v, &surface);
    if(!success) {
        //print_vector(&v);
        //printf("%lf\n", eval_ext(bargs->q, &v));
        //print_vector(&surface);
        //printf("%lf\n", eval_ext(bargs->q, &surface));
        exit(1);
    }
    bargs->func(bargs->s, bargs->q, &surface); 
}

void multi_thread_benchmark(int64_t radius, int num_threads) {
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -radius * radius};
    pthread_t *threads = 
        (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    builder_args *args = 
        (builder_args *)malloc(num_threads * sizeof(builder_args));

    struct timespec start, end;
    int64_t elapsed;

    uint64_t i, j;
    int trials = 1;
    subspace *s;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i = 0; i < trials; i++) {
        s = subspace_init(-radius - 1, -radius - 1, 0, radius + 2,
            radius + 2, 1);
        uint64_t vol = volume(s);
        for (j = 0; j < num_threads; j++) {
            args[j].s = s;
            args[j].q = &q;
            args[j].func = breadth_first_fill;
            args[j].index = j * vol / num_threads;
            pthread_create(&threads[j], NULL, builder_thread, &args[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        //display_subspace(s);
        subspace_free(s);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = 1000000000 * (uint64_t)(end.tv_sec - start.tv_sec) +
            (uint64_t)end.tv_nsec - (uint64_t)start.tv_nsec;
    
    printf("Elapsed time for breadth first search is: %lld nanoseconds\n",
            elapsed);

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i = 0; i < trials; i++) {
        s = subspace_init(-radius - 1, -radius - 1, 0, radius + 2,
            radius + 2, 1);
        uint64_t vol = volume(s);
        for (j = 0; j < num_threads; j++) {
            args[j].s = s;
            args[j].q = &q;
            args[j].func = depth_first_fill;
            args[j].index = j * vol / num_threads;
            pthread_create(&threads[j], NULL, builder_thread, &args[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        //display_subspace(s);
        subspace_free(s);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = 1000000000 * (uint64_t)(end.tv_sec - start.tv_sec) +
            (uint64_t)end.tv_nsec - (uint64_t)start.tv_nsec;
    
    printf("Elapsed time for depth first search is: %lld nanoseconds\n",
            elapsed);

    free(threads);
    free(args);
}

void single_thread_benchmark(int64_t radius) {
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -radius * radius};
    int64_t i, trials = 100;
    struct timespec start, end;
    subspace *s;
    vector v = {0, 0, 0};
    int64_t elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i = 0; i < trials; i++) {
        s = subspace_init(-radius - 1, -radius - 1, 0, radius + 2,
            radius + 2, 1);
        breadth_first_fill(s, &q, &v);
        subspace_free(s);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = 1000000000 * (uint64_t)(end.tv_sec - start.tv_sec) +
            (uint64_t)end.tv_nsec - (uint64_t)start.tv_nsec;
    
    printf("Elapsed time for breadth first search is: %lld nanoseconds\n",
            elapsed);

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i = 0; i < trials; i++) {
        s = subspace_init(-radius - 1, -radius - 1, 0, radius + 2,
            radius + 2, 1);
        depth_first_fill(s, &q, &v);
        subspace_free(s);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = 1000000000 * (uint64_t)(end.tv_sec - start.tv_sec) +
            (uint64_t)end.tv_nsec - (uint64_t)start.tv_nsec;

    printf("Elapsed time for depth first search is: %lld nanoseconds\n",
            elapsed);

}

void display_subspace(subspace *s) {
    int64_t i, j, k;
    for (k = s->z_min; k < s->z_max; k++) {
        for (j = s->y_max - 1; j >= s->y_min; j--) {
            for (i = s->x_min; i < s->x_max; i++) {
                printf("%d ", s->points[_index(s, i, j, k)].surface);
            }
            putchar('\n');
        }
    }
}
