#ifndef QUADRIC_H
#define QUADRIC_H
#include <stdint.h>
#include <semaphore.h>

/* Quadratic surfaces are also called quadrics, and there are 17 
 * standard-form types. A quadratic surface intersects every plane in a 
 * (proper or degenerate) conic section. In addition, the cone consisting 
 * of all tangents from a fixed point to a quadratic surface cuts every plane 
 * in a conic section, and the points of contact of this cone with the 
 * surface form a conic section (Hilbert and Cohn-Vossen 1999, p. 12). */

/* ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0. */

#define _index(s, x, y, z) (\
    ((x) - (s)->x_min) * ((s)->y_max - (s)->y_min) * ((s)->z_max - (s)->z_min) + \
    ((y) - (s)->y_min) * ((s)->z_max - (s)->z_min) + \
    ((z) - (s)->z_min)\
)

#define frozen_point() ( )

#define _x(s, index) ((index) / ((s)->y_max - (s)->y_min) / \
        ((s)->z_max - (s)->z_min) + (s)->x_min)
#define _y(s, index) ((index) % (((s)->y_max - (s)->y_min) * \
        ((s)->z_max - (s)->z_min)) / ((s)->z_max - (s)->z_min) + (s)->y_min)
#define _z(s, index) ((index) % (((s)->y_max - (s)->y_min) * \
        ((s)->z_max - (s)->z_min)) % ((s)->z_max - (s)->z_min) + (s)->z_min)

#define touch(s) (sem_post(&(s)->points_plotted))

#define volume(s) ((size_t)(((s)->x_max - (s)->x_min) * \
        ((s)->y_max - (s)->y_min) * \
        ((s)->z_max - (s)->z_min)))


typedef struct _quadric {
    double a, b, c, d, e, f, g, h, i, j;
} quadric;

typedef struct _vector {
    double x, y, z;
} vector;

typedef struct _point {
    int64_t x, y, z;
    sem_t sema;
    uint8_t plotted;
} point;

typedef struct _subspace {
    int64_t x_min, y_min, z_min, x_max, y_max, z_max; 
    sem_t points_plotted;
    point *points;
} subspace;

typedef struct _frozen_subspace {
    int64_t x_min, y_min, z_min, x_max, y_max, z_max; 
    /* Bit field for each point in the bounding volume. 1 if it is plotted,
     * 0 otherwise */
    uint8_t *points;
} frozen_subspace;

typedef struct _node {
    void *data;
    struct _node *prev;
    struct _node *next;
} node;

typedef struct _list {
    node begin;
    node end;
} list;

subspace *subspace_init(int64_t, int64_t, int64_t, 
        int64_t, int64_t, int64_t);
void subspace_free(subspace *);
void frozen_subspace_free(frozen_subspace *);
double eval_int(const quadric *, const vector *);
double eval_ext(const quadric *, const vector *);
int is_surface(const quadric *, const vector *);
void print_vector(const vector *v);
void print_subspace(const subspace *s);
int find_surface(quadric *q, const vector *, vector *);
void depth_first_surface(subspace *, const quadric *, const vector *);
void breadth_first_surface(subspace *, const quadric *, const vector *);
void depth_first_fill(subspace *, const quadric *, const vector *);
void breadth_first_fill(subspace *, const quadric *, const vector *);
list *new_list();
void *pop(list *);
void *peek(list *);
void *push_back(list *, void *);
void list_destroy(list *);
void print_list(list *, void (*)(void *));
int empty(list *);

double eval_ext(const quadric *, const vector *);
#endif
