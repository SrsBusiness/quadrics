#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <semaphore.h>
#include <pthread.h>
#include <errno.h>
#include "minicurses.h"

/* Quadratic surfaces are also called quadrics, and there are 17 
 * standard-form types. A quadratic surface intersects every plane in a 
 * (proper or degenerate) conic section. In addition, the cone consisting 
 * of all tangents from a fixed point to a quadratic surface cuts every plane 
 * in a conic section, and the points of contact of this cone with the 
 * surface form a conic section (Hilbert and Cohn-Vossen 1999, p. 12). */

/* ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0. */

#define index(s, x, y, z) (\
    ((x) - (s)->x_min) * ((s)->y_max - (s)->y_min) * ((s)->z_max - (s)->z_min) + \
    ((y) - (s)->y_min) * ((s)->z_max - (s)->z_min) + \
    ((z) - (s)->z_min)\
)

#define x(s, index) ((index) / ((s)->y_max - (s)->y_min) / \
        ((s)->z_max - (s)->z_min) + (s)->x_min)
#define y(s, index) ((index) % (((s)->y_max - (s)->y_min) * \
        ((s)->z_max - (s)->z_min)) / ((s)->z_max - (s)->z_min) + (s)->y_min)
#define z(s, index) ((index) % (((s)->y_max - (s)->y_min) * \
        ((s)->z_max - (s)->z_min)) % ((s)->z_max - (s)->z_min) + (s)->z_min)


typedef struct _quadric {
    double a, b, c, d, e, f, g, h, i, j;
} quadric;

typedef struct _vector {
    int64_t x, y, z;
} vector;

typedef struct _point {
    int64_t x, y, z;
    sem_t sema;
    int surface;
} point;

typedef struct _subspace {
    int64_t x_min, y_min, z_min, x_max, y_max, z_max; 
    point *points;
} subspace;

subspace *subspace_init(uint64_t, uint64_t, uint64_t, 
        uint64_t, uint64_t, uint64_t);
void subspace_free(subspace *);
double eval_int(quadric *, vector *);
double eval_ext(quadric *, vector *);
int is_surface(quadric *, vector *);
void print_vector(vector *v);
void print_subspace(subspace *s);
