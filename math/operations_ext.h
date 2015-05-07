#pragma once

#include "types.h"
#include "operations.h"
#include <math.h>

namespace math {
    inline void operator /=(vec3& a, point_t b) {
        a.x /= b;
        a.y /= b;
        a.z /= b;
    }

    inline void normalize(vec3& vec) {
        auto norm = sqrtf(dot(vec, vec));
        if (norm > 0) {
            vec /= norm;
        }
    }
}