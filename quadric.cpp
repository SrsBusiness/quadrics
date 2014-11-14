#include "quadric.h"
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <boost/circular_buffer.hpp>

#define EPSILON 0.5
/* Quadratic surfaces are also called quadrics, and there are 17 
 * standard-form types. A quadratic surface intersects every plane in a 
 * (proper or degenerate) conic section. In addition, the cone consisting 
 * of all tangents from a fixed point to a quadratic surface cuts every plane 
 * in a conic section, and the points of contact of this cone with the 
 * surface form a conic section (Hilbert and Cohn-Vossen 1999, p. 12). */

/* ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0. */

double (*eval[])(quadric *, vector *) = {eval_ext, eval_int};

double eval_int(quadric *q, vector *v) {
    double result = q->a * v->x * v->x + 
        q->b * v->y * v->y + 
        q->c * v->z * v->z + 
        q->d * v->y * v->z + 
        q->e * v->x * v->z + 
        q->f * v->x * v->y+ 
        q->g * v->x + 
        q->h * v->y + 
        q->i * v->z + 
        q->j;
    // bias of 0 should go to negative
    return result == 0.0 ? -result : result;
}

double eval_ext(quadric *q, vector *v) {
    double result = q->a * v->x * v->x + 
        q->b * v->y * v->y + 
        q->c * v->z * v->z + 
        q->d * v->y * v->z + 
        q->e * v->x * v->z + 
        q->f * v->x * v->y+ 
        q->g * v->x + 
        q->h * v->y + 
        q->i * v->z + 
        q->j;
    // bias of 0 should go to positive
    return result;
}

/* if eval() of all neighboring points are of all the same sign, then
 * v cannot be on the surface. Otherwise, it may be*/
int is_surface(quadric *q, vector *v) {
    double val = eval_ext(q, v);
    /* short circuit */
    if (val == 0.0)
        return 1;
    uint64_t inside = (0x8000000000000000 & *(uint64_t *)&val) >> 63;
    double (*eval_func)(quadric *, vector *) = eval[inside];
    vector tmp;
    int i = 1;
    uint64_t sign1, sign2;
    for (; i <= 4; i <<= 1) {
        tmp.x = v->x + EPSILON * (i & 0x1);
        tmp.y = v->y + EPSILON * ((i & 0x2) >> 1);
        tmp.z = v->z + EPSILON * ((i & 0x4) >> 2);
        val = eval_func(q, &tmp);
        sign1 = 0x8000000000000000 & *(uint64_t *)&val;
        tmp.x = v->x - EPSILON * (i & 0x1);
        tmp.y = v->y - EPSILON * ((i & 0x2) >> 1);
        tmp.z = v->z - EPSILON * ((i & 0x4) >> 2);
        val = eval_func(q, &tmp);
        sign2 = 0x8000000000000000 & *(uint64_t *)&val;
        if(sign1 != sign2)
            return 1;
    }
    return 0;
}

subspace *subspace_init(uint64_t x_min, uint64_t y_min, 
        uint64_t z_min, uint64_t x_max, uint64_t y_max, uint64_t z_max) {
    subspace *s = (subspace *)malloc(sizeof(subspace));
    s->x_min = x_min;
    s->y_min = y_min;
    s->z_min = z_min;
    s->x_max = x_max;
    s->y_max = y_max;
    s->z_max = z_max;
    size_t space_size = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    s->points = (point *)malloc(space_size * sizeof(point));
    size_t i;
    for (i = 0; i < space_size; i++) {
        // initialize x, y, and z
        s->points[i].x = _x(s, i);
        s->points[i].y = _y(s, i);
        s->points[i].z = _z(s, i);
        if(sem_init(&s->points[i].sema, 0, 1) == -1) {
            free(s->points);
            free(s);
            return NULL;
        }
    }
    return s;
}

size_t volume(subspace *s) {
    return (s->x_max - s->x_min) * 
        (s->y_max - s->y_min) * 
        (s->z_max - s->z_min);
}

void subspace_free(subspace *s) {
    free(s->points);
    free(s);
}

void depth_first_fill(subspace *s, quadric *q, vector *v) {
    /* if out of bounding volume */
    if (v->x < s->x_min || v->x >= s->x_max ||
            v->y < s->y_min || v->y >= s->y_max ||
            v->z < s->z_min || v->z >= s->z_max)
        return;

    uint64_t index = _index(s, v->x, v->y, v->z);
    /* if already visited */
    if(sem_trywait(&s->points[index].sema) == -1 &&
            errno == EAGAIN)
        return;
    
    s->points[index].surface = is_surface(q, v);
    int i;
    vector tmp;
    for (i = 0; i <= 4; i <<=1) {
        tmp.x = v->x + (i & 0x1);
        tmp.y = v->y + ((i & 0x2) >> 1);
        tmp.z = v->z + ((i & 0x4) >> 2);
        depth_first_fill(s, q, &tmp);
        tmp.x = v->x - (i & 0x1);
        tmp.y = v->y - ((i & 0x2) >> 1);
        tmp.z = v->z - ((i & 0x4) >> 2);
        depth_first_fill(s, q, &tmp);
    }
}

void breadth_first_fill(subspace *s, quadric *q, vector *v) {
    boost::circular_buffer<vector *> queue = 
        boost::circular_buffer<vector *>(volume(s));

    vector *tmp, *current = (vector *)malloc(sizeof(vector));
    current->x = v->x; current->y = v->y; current->z = v->z; 
    queue.push_back(current);

    int i;
    while (!queue.empty()) {
        current = queue.front();
        queue.pop_front();
    
        if (current->x < s->x_min || current->x >= s->x_max ||
            current->y < s->y_min || current->y >= s->y_max ||
            current->z < s->z_min || current->z >= s->z_max)
            continue;

        uint64_t index = _index(s, current->x, current->y, current->z);
        if(sem_trywait(&s->points[index].sema) == -1 &&
            errno == EAGAIN)
            continue;
        s->points[index].surface = is_surface(q, current);
        
        for (i = 1; i <= 4; i <<=1) {
            tmp = (vector *)malloc(sizeof(vector));
            tmp->x = current->x + (i & 0x1);
            tmp->y = current->y + ((i & 0x2) >> 1);
            tmp->z = current->z + ((i & 0x4) >> 2);
            queue.push_back(tmp);
            tmp = (vector *)malloc(sizeof(vector));
            tmp->x = current->x - (i & 0x1);
            tmp->y = current->y - ((i & 0x2) >> 1);
            tmp->z = current->z - ((i & 0x4) >> 2);
            queue.push_back(tmp);
        }
        free(current);
    }
    return;
}

void print_subspace(subspace *s) {
    printf("x: %lld to %lld, y: %lld to %lld, z: %lld to %lld\n",
            s->x_min, s->x_max, s->y_min, s->y_max, s->z_min, s->z_max);
}

void print_vector(vector *v) {
    printf("{%lld, %lld, %lld}\n", v->x, v->y, v->z);
}
