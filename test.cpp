#include "quadric.h"
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "minicurses.h"
#include <pthread.h>
#include <assert.h>
#include <math.h>
#include "serialize.h"


void display_subspace(subspace *);
void display_frozen_subspace(frozen_subspace *);

void single_thread_benchmark(int64_t);
void multi_thread_benchmark(int64_t, int);
void find_surface_test(int64_t);
void herp_test();
void serialize_test(int64_t);

int main(int argc, char **argv) {
    //multi_thread_benchmark(19, 32);
    //single_thread_benchmark(19);
    find_surface_test(18);
    //herp_test();
    //serialize_test(64);
}

void serialize_test(int64_t radius) {
    subspace *s = subspace_init(-radius - 1, -radius - 1, 
            -radius - 1, radius + 2, radius + 2, radius + 2);
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -radius * radius};
    vector surface, v = {0, 0, 0};
    find_surface(&q, &v, &surface);
    breadth_first_surface(s, &q, &surface);
    int points;
    sem_getvalue(&s->points_plotted, &points);
    size_t len = volume(s);
    uint8_t *buf = (uint8_t *)calloc((len + 7) / 8, sizeof(uint8_t));
    subspace_dump(s, buf, 0, len);
    FILE *f = fopen("circle64_expected", "w");
    fwrite(buf, sizeof(char), (len + 7) / 8, f);
    fclose(f);
        subspace_serialize(s, "circle64.lzma");
    frozen_subspace *fs = frozen_subspace_deserialize("circle64.lzma", 
            0, 0, 0); 
    int64_t i, j;
    printf("Subspace:\n");
    for (i = 0; i < volume(s); i+= 64) {
        for (j = (i + 63) < (volume(s) - 1) ? 63 : (volume(s) - 1) % 64; j >= 0; j--) {
            putchar('0' + s->points[i + j].plotted);
        }
        putchar('\n');
    }
    printf("\nBuffer:\n");
    for (i = 0; i < (len  + 7) / 8; i += sizeof(uint64_t)) {
        printf("%llx\n", *(uint64_t *)(buf + i));
    }

    frozen_subspace *s_frozen = freeze_subspace(s);
    FILE *fsub = fopen("frozen_subspace", "w");
    fwrite(fs->points, 1, ((fs->x_max - fs->x_min) * (fs->y_max - fs->y_min) * 
            (fs->z_max - fs->z_min) + 7) / 8, fsub);
    fclose(fsub);
    frozen_subspace_serialize(fs, "frozen_subspace_serialize.lzma");
    frozen_subspace_serialize(fs, "s_frozen_serialize.lzma");
    frozen_subspace_free(fs);
    subspace_free(s);
    frozen_subspace_free(s_frozen);
    free(buf);
}

void herp_test() {
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -19 * 19};
    vector v = {-10, 16, 2};
    is_surface(&q, &v);
}

void find_surface_test(int64_t radius) {
    subspace *s = subspace_init(-radius - 1, -radius - 1, -radius - 1, radius + 2,
            radius + 2, radius + 2);
    quadric q = {1, 1, 1, 0, 0, 0, 0, 0, 0, -radius * radius};
    vector v = {0, 0, 0};
    vector surface;
    assert(find_surface(&q, &v, &surface)); 
    depth_first_fill(s, &q, &surface);
    subspace_serialize(s, "hello");
    //frozen_subspace *f = freeze_subspace(s);
    frozen_subspace *f = frozen_subspace_deserialize("hello", s->x_min, s->y_min, s->z_min);
    display_frozen_subspace(f);
}

void print_elem(void *elem) {
    printf("%ld\n", (int64_t)elem);
}

typedef struct _builder_args {
    subspace *s;
    quadric *q;
    int64_t index;
    void (*func)(subspace *, const quadric *, const vector *);
} builder_args;

void *builder_thread(void *args) {
    builder_args *bargs = (builder_args *)args;
    vector v, surface;
    v.x = _x(bargs->s, bargs->index);
    v.y = _y(bargs->s, bargs->index);
    v.z = _z(bargs->s, bargs->index);
    if(!find_surface(bargs->q, &v, &surface)) {
        print_vector(&surface);
        exit(1);
    }
    bargs->func(bargs->s, bargs->q, &surface); 
}

void multi_thread_benchmark(int64_t radius, int num_threads) {
    quadric q = {1, -1, 0, 0, 0, 0, 0, 0, -1, 0};
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
        s = subspace_init(-radius - 1, -radius - 1, -radius - 1, radius + 2,
            radius + 2, radius + 2);
        uint64_t vol = volume(s);
        for (j = 0; j < num_threads; j++) {
            args[j].s = s;
            args[j].q = &q;
            args[j].func = breadth_first_surface;
            args[j].index = j * vol / num_threads;
            pthread_create(&threads[j], NULL, builder_thread, &args[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        display_subspace(s);
        int points_plotted;
        sem_getvalue(&s->points_plotted, &points_plotted);
        printf("%d points plotted\n", points_plotted);
        subspace_free(s);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = 1000000000 * (uint64_t)(end.tv_sec - start.tv_sec) +
            (uint64_t)end.tv_nsec - (uint64_t)start.tv_nsec;
    
    printf("Elapsed time for breadth first search is: %lld nanoseconds\n",
            elapsed);

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i = 0; i < trials; i++) {
        s = subspace_init(-radius - 1, -radius - 1, -radius - 1, radius + 2,
            radius + 2, radius + 2);
        uint64_t vol = volume(s);
        for (j = 0; j < num_threads; j++) {
            args[j].s = s;
            args[j].q = &q;
            args[j].func = depth_first_surface;
            args[j].index = j * vol / num_threads;
            pthread_create(&threads[j], NULL, builder_thread, &args[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        display_subspace(s);
        int points_plotted;
        sem_getvalue(&s->points_plotted, &points_plotted);
        printf("%d points plotted\n", points_plotted);
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
        breadth_first_surface(s, &q, &v);
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
        depth_first_surface(s, &q, &v);
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
        clear_all();
        for (j = s->y_max - 1; j >= s->y_min; j--) {
            for (i = s->x_min; i < s->x_max; i++) {
                printf("%d ", s->points[_index(s, i, j, k)].plotted);
            }
            putchar('\n');
        }
        putchar('\n');
        sleep(1);
    }
}

void display_frozen_subspace(frozen_subspace *s) {
    printf("%d\n", s->z_min);
    int64_t i, j, k;
    for (k = s->z_min; k < s->z_max; k++) {
        clear_all();
        for (j = s->y_max - 1; j >= s->y_min; j--) {
            for (i = s->x_min; i < s->x_max; i++) {
                printf("%d ", frozen_point(s, _index(s, i, j, k)));
            }
            putchar('\n');
        }
        putchar('\n');
        sleep(1);
    }
}
