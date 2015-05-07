#pragma once

#include "types.h"

namespace math {

    struct plane_t {
        vec3 origin;
        vec3 normal;
    };

    point_t intr_ray_plane(plane_t& plane, ray_t& ray, vec3& intr);
}