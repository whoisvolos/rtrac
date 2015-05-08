#pragma once

#include "types.h"
#include "operations.h"
#include <math.h>

#define EPSILON   0.00000001

namespace math {
    struct matrix3x3 {
        point_t p[3][3];
    };

    void dump(vec3& vec);
    void dump(matrix3x3& mat);

    inline void operator /=(vec3& a, point_t b) {
        a.x /= b;
        a.y /= b;
        a.z /= b;
    }

    inline point_t norm(vec3& vec) {
        return sqrtf(dot(vec, vec));
    }

    inline void normalize(vec3& vec) {
        point_t n = norm(vec);
        if (n != 0) {
            vec /= sqrt(n);
        }
    }

    inline matrix3x3 make_mat3x3(point_t a00, point_t a01, point_t a02, point_t a10, point_t a11, point_t a12, point_t a20, point_t a21, point_t a22)
    {
        return { { { a00, a01, a02 }, { a10, a11, a12 }, { a20, a21, a22 } } };
    }

    matrix3x3 IDENTITY_33 = make_mat3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    matrix3x3 ZERO_33 = make_mat3x3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    vec3 X = make_vec3(1, 0, 0);
    vec3 Y = make_vec3(0, 1, 0);
    vec3 Z = make_vec3(0, 0, 1);

    inline vec3 operator * (matrix3x3& mat, vec3& vec) {
        return make_vec3(
            (mat.p[0][0] * vec.x) + (mat.p[1][0] * vec.y) + (mat.p[2][0] * vec.z),
            (mat.p[0][1] * vec.x) + (mat.p[1][1] * vec.y) + (mat.p[2][1] * vec.z),
            (mat.p[0][2] * vec.x) + (mat.p[1][2] * vec.y) + (mat.p[2][2] * vec.z)
        );
    }

    inline matrix3x3 operator * (matrix3x3& a, matrix3x3& b) {
        return make_mat3x3(
            a.p[0][0] * b.p[0][0] + a.p[0][1] * b.p[1][0] + a.p[0][2] * b.p[2][0],
            a.p[1][0] * b.p[0][0] + a.p[1][1] * b.p[1][0] + a.p[1][2] * b.p[2][0],
            a.p[2][0] * b.p[0][0] + a.p[2][1] * b.p[1][0] + a.p[2][2] * b.p[2][0],
            a.p[0][0] * b.p[0][1] + a.p[0][1] * b.p[1][1] + a.p[0][2] * b.p[2][1],
            a.p[1][0] * b.p[0][1] + a.p[1][1] * b.p[1][1] + a.p[1][2] * b.p[2][1],
            a.p[2][0] * b.p[0][1] + a.p[2][1] * b.p[1][1] + a.p[2][2] * b.p[2][1],
            a.p[0][0] * b.p[0][2] + a.p[0][1] * b.p[1][2] + a.p[0][2] * b.p[2][2],
            a.p[1][0] * b.p[0][2] + a.p[1][1] * b.p[1][2] + a.p[1][2] * b.p[2][2],
            a.p[2][0] * b.p[0][2] + a.p[2][1] * b.p[1][2] + a.p[2][2] * b.p[2][2]
        );
    }

    inline matrix3x3 ssc_mul (matrix3x3& a, matrix3x3& b) {
        return make_mat3x3(
            a.p[0][1] * b.p[1][0] + a.p[0][2] * b.p[2][0],
            a.p[1][2] * b.p[2][0],
            a.p[2][1] * b.p[1][0],
            a.p[0][1] * b.p[1][1] + a.p[0][2] * b.p[2][1],
            a.p[1][0] * b.p[0][1] + a.p[1][2] * b.p[2][1],
            a.p[2][0] * b.p[0][1],
            a.p[0][1] * b.p[1][2],
            a.p[1][0] * b.p[0][2],
            a.p[2][0] * b.p[0][2] + a.p[2][1] * b.p[1][2]
        );
    }

    inline matrix3x3 operator + (matrix3x3& a, matrix3x3& b) {
        return make_mat3x3(
            a.p[0][0] + b.p[0][0], a.p[1][0] + b.p[1][0], a.p[2][0] + b.p[2][0],
            a.p[0][1] + b.p[0][1], a.p[1][1] + b.p[1][1], a.p[2][1] + b.p[2][1],
            a.p[0][2] + b.p[0][2], a.p[1][2] + b.p[1][2], a.p[2][2] + b.p[2][2]
        );
    }

    inline matrix3x3 operator * (matrix3x3& a, point_t b) {
        return make_mat3x3(
            a.p[0][0] * b, a.p[1][0] * b, a.p[2][0] * b,
            a.p[0][1] * b, a.p[1][1] * b, a.p[2][1] * b,
            a.p[0][2] * b, a.p[1][2] * b, a.p[2][2] * b
        );
    }

    inline matrix3x3 rotate_towards(vec3 subject, vec3 to) {
        // TODO: Normalize?
        point_t c = dot(subject, to);
        // Parallel check
        if (fabsf(c - 1) >= EPSILON) {
            vec3 v = cross(subject, to);
            point_t s2 = dot(v, v);
            matrix3x3 ssc = make_mat3x3(0, v.z, -v.y, -v.z, 0, v.x, v.y, -v.x, 0);
            matrix3x3 ssc2 = ssc_mul(ssc, ssc);
            ssc2 = ssc2 * ((1 - c) / s2);
            matrix3x3 rot = IDENTITY_33 + ssc;
            rot = rot + ssc2;
            return rot;
        } else {
            return IDENTITY_33;
        }
    }

    void dump(vec3& vec) {
        printf("%f\t%f\t%f\n", vec.x, vec.y, vec.z);
    }

    void dump(matrix3x3& mat) {
        printf("%f\t%f\t%f\n", mat.p[0][0], mat.p[1][0], mat.p[2][0]);
        printf("%f\t%f\t%f\n", mat.p[0][1], mat.p[1][1], mat.p[2][1]);
        printf("%f\t%f\t%f\n", mat.p[0][2], mat.p[1][2], mat.p[2][2]);
    }
}