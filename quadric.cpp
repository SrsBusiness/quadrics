#include "quadric.h"
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <list>
#include <float.h>
#include <math.h>
//#include <boost/circular_buffer.hpp>

#define EPSILON 0.50
/* Quadratic surfaces are also called quadrics, and there are 17 
 * standard-form types. A quadratic surface intersects every plane in a 
 * (proper or degenerate) conic section. In addition, the cone consisting 
 * of all tangents from a fixed point to a quadratic surface cuts every plane 
 * in a conic section, and the points of contact of this cone with the 
 * surface form a conic section (Hilbert and Cohn-Vossen 1999, p. 12). */

/* ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0. */

double (*eval_funcs[])(const quadric *, const vector *) = {eval_ext, eval_int};
double (*eval)(const quadric *, const vector *) = eval_ext;

double eval_int(const quadric *q, const vector *v) {
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

double eval_ext(const quadric *q, const vector *v) {
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
int is_surface(const quadric *q, const vector *v) {
    double val = eval_ext(q, v);
    if (val == 0.0)
        return 1;
    vector tmp;
    int i = 1;
    uint64_t sign1, sign2;
    for (; i <= 4; i <<= 1) {
        tmp.x = v->x + EPSILON * (i & 0x1);
        tmp.y = v->y + EPSILON * ((i & 0x2) >> 1);
        tmp.z = v->z + EPSILON * ((i & 0x4) >> 2);
        val = eval_ext(q, &tmp);
        if (val == 0.0) 
            return 0;
        sign1 = val > 0.0;
        tmp.x = v->x - EPSILON * (i & 0x1);
        tmp.y = v->y - EPSILON * ((i & 0x2) >> 1);
        tmp.z = v->z - EPSILON * ((i & 0x4) >> 2);
        val = eval_ext(q, &tmp);
        sign2 = val > 0.0;
        if (val == 0.0)
            return 0;
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
    s->points = (point *)calloc(space_size, sizeof(point));
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

/* Given a starting point v, traces a path until it arrives at a surface
 * point. Returns 1 if successful, 0 otherwise. Surface is modified with the
 * coordinates of the surface point. If 0 is returned, the new coordinates in 
 * surface are undefined
 */
int find_surface(quadric *q, const vector *v, vector *surface) {
    vector tmp, closest;
    double dist, shortest_dist;
    int i;
    shortest_dist = DBL_MAX;
    int progress = 1;
    *surface = *v;
    while(!is_surface(q, surface) && progress) {
        progress = 0;
        for (i = 1; i <= 4; i <<=1) {
            tmp.x = surface ->x + (i & 0x1);
            tmp.y = surface ->y + ((i & 0x2) >> 1);
            tmp.z = surface ->z + ((i & 0x4) >> 2);
            if ((dist = fabs(eval_ext(q, &tmp))) < shortest_dist) {
                shortest_dist = dist;
                closest = tmp;
                progress = 1;
            }
            tmp.x = surface->x - (i & 0x1);
            tmp.y = surface->y - ((i & 0x2) >> 1);
            tmp.z = surface->z - ((i & 0x4) >> 2);
            if ((dist = fabs(eval_ext(q, &tmp))) < shortest_dist) {
                shortest_dist = dist;
                closest = tmp;
                progress = 1;
            }
        }
        *surface = closest;
    }
    return progress;
}

/* precondition: v is surface point
 * Depth first trace of all surface points
 * */
void depth_first_fill(subspace *s, const quadric *q, const vector *v) {
    /* if out of bounding volume */
    if (v->x < s->x_min || v->x >= s->x_max ||
            v->y < s->y_min || v->y >= s->y_max ||
            v->z < s->z_min || v->z >= s->z_max)
        return;

    uint64_t index = _index(s, v->x, v->y, v->z);
    /* if already visited */
    if (sem_trywait(&s->points[index].sema) == -1 &&
            errno == EAGAIN)
        return;
   
    /* if not a surface point */
    if (!(s->points[index].surface = is_surface(q, v)))
        return;

    vector tmp;
    int i, j, k;
    for (i = -1; i <= 1; i++) {
        for (j = -1; j <= 1; j++) {
            for (k = -1; k <= 1; k++) {
                if (!i && !j && !k)
                    continue;
                tmp.x = v->x + i;
                tmp.y = v->y + j;
                tmp.z = v->z + k;
                depth_first_fill(s, q, &tmp);
            }
        }
    }
}

void print_func(void *data) {
    printf("%p\n", data); 
}

void breadth_first_fill(subspace *s, const quadric *q, const vector *v) {
    std::list<vector *> queue = std::list<vector *>();
    vector *tmp, *current = (vector *)malloc(sizeof(vector));
    current->x = v->x; current->y = v->y; current->z = v->z; 

    queue.push_back(current);
    int i, j, k;
    uint64_t surface_points = 0;
    while (!queue.empty()) {
        current = queue.front();
        queue.pop_front();
        uint64_t index;
        if (current->x < s->x_min || current->x >= s->x_max ||
            current->y < s->y_min || current->y >= s->y_max ||
            current->z < s->z_min || current->z >= s->z_max)
            goto cleanup;

        index = _index(s, current->x, current->y, current->z);
        if(sem_trywait(&s->points[index].sema) == -1 &&
            errno == EAGAIN)
            goto cleanup;

        /* If point is not on surface, do not visit */
        if (!(s->points[index].surface = is_surface(q, current)))
            goto cleanup;
        //else {
        //    surface_points++;
        //}
        for (i = -1; i <= 1; i++) {
            for (j = -1; j <= 1; j++) {
                for (k = -1; k <= 1; k++) {
                    if (!i && !j && !k)
                        continue;
                    tmp = (vector *)malloc(sizeof(vector));
                    tmp->x = current->x + i;
                    tmp->y = current->y + j;
                    tmp->z = current->z + k;
                    queue.push_back(tmp);
                }
            }
        }
        cleanup:
        free(current);
    }
    printf("%llu points plotted\n", surface_points);
    return;
}

void print_subspace(const subspace *s) {
    printf("x: %lld to %lld, y: %lld to %lld, z: %lld to %lld\n",
            s->x_min, s->x_max, s->y_min, s->y_max, s->z_min, s->z_max);
}

void print_vector(const vector *v) {
    printf("{%lf, %lf, %lf}\n", v->x, v->y, v->z);
}

