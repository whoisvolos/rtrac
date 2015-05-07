#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "funcs.hpp"
#include "time.hpp"

#include <random>
#include "math/types.h"
#include "math/plane.h"

#include <assert.h>

using namespace math;

template <class GEN>
vec3 nrand_hemisphere(GEN &gen_rad, GEN &gen_phi) {
    float r = gen_rad();
    float rad = sqrtf(r);
    float phi = gen_phi() * 2 * M_PI;
    return { rad * cosf(phi), rad * sinf(phi), sqrtf(1 - r) };
}

template <class GEN>
vec3 nrand_hemisphere_spherical(GEN &gen_rad, GEN &gen_phi) {
    float theta = acosf(sqrtf(gen_rad()));
    float phi = gen_phi() * 2 * M_PI;
    return vec3 { cosf(phi) * sinf(theta), sinf(phi) * sinf(theta), cosf(theta) };
};

/**
 * Parallel plates
 * param a - x dimension
 * param b - y dimension
 * param c - distance between plates
 */
point_t real_result(point_t a, point_t b, point_t c) {
    auto x = a/c, y = b/c, xq = x * x, yq = y * y;
    auto result = 2 / M_PI / x / y * (logf(sqrtf((1 + xq) * (1 + yq) / (1 + xq + yq))) + x * sqrtf(1 + yq) * atanf(x / sqrtf(1 + yq)) + y * sqrtf(1 + xq) * atanf(y / sqrtf(1 + xq)) - x * atanf(x) - y * atanf(y));
    return result;
}

/**
 * Perp plates
 * param a - 1st plate x width
 * param b - Both plates y width
 * param c - 2nd plane z height
 */
double real_result_90(double a, double b, double c) {
    double H = c / b;
    double W = a / b;
    double WH2 = H * H + W * W;

    double ln_part = 0.25 * log( ((1 + W * W)*(1 + H * H) / (1 + WH2)) * pow((W * W * (1 + WH2) / (1 + W * W) / WH2), W * W) * pow((H * H * (1 + WH2) / (1 + H * H) / WH2), H * H) );
    double result = (1 / M_PI / W) * (W * atan(1 / W) + H * atan(1 / H) - sqrt(WH2) * atan(1 / sqrt(WH2)) + ln_part);

    return result;
}


int main(int argc, char **argv) {
    const int TOTAl_RAYS = 1000;
    float c = 1, a = 1, b = 1;

    plane_t A1 = { { 0, 0, 0 }, { 0, 0, 1 } };
    plane_t A2 = { { 0, 0, c }, { 1, 0, 0 } };

    std::mt19937 eng1, eng2(1000);
    std::uniform_real_distribution<float> x_distr(0, 1);

    auto rgen = std::bind(x_distr, eng1);
    auto pgen = std::bind(x_distr, eng2);
    unsigned long long ok = 0;

    struct timespec t2, t3;
    clock_gettime(CLOCK_MONOTONIC,  &t2);

    for (int i = 0; i < TOTAl_RAYS; ++i) {
        // Generate random position
        vec3 p0 = { a * rgen(), b * pgen(), 0 };

        // Generate cosine-weighted random ray
        auto u = nrand_hemisphere(rgen, pgen);
        ray_t ray = { p0, u };

        // Calculate intersection with A2 plane
        vec3 intr;
        float dist = intr_ray_plane(A2, ray, intr);
        if (dist > 0 &&
            intr.z >= 0 && intr.z < c &&
            intr.y >= 0 && intr.y < b) {
            ++ok;
        }

    };
    auto result = (float)ok/(float)(TOTAl_RAYS);

    clock_gettime(CLOCK_MONOTONIC,  &t3);
    auto dt1 = (t3.tv_sec - t2.tv_sec) + (float) (t3.tv_nsec - t2.tv_nsec) * 1e-9;

    auto real = real_result_90(a, b, c);
    printf("Result: %f, real: %f, rel: %f, rays: %i\n", result, real, result / real, TOTAl_RAYS);
    printf("%.3f ms\n", dt1 * 1000);

    return 0;
}
