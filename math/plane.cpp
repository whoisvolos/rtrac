#include "plane.h"
#include "operations.h"
#include <math.h>

#define EPSILON   0.00000001

namespace math {
    point_t intr_ray_plane(plane_t& plane, ray_t& ray, vec3& intr) {
        point_t denom = dot(plane.normal, ray.direction);
        if (fabs(denom) < EPSILON) {
            return -1;
        }
        point_t s = dot(plane.normal, (plane.origin - ray.origin)) / denom;
        intr = ray.origin + ray.direction * s;
        return s;
    }
}