#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "funcs.hpp"
#include "time.hpp"

#include <random>
#include "linalg.hpp"

#include <assert.h>

using namespace linalg;

//  globals
int numnodes, myid, mpi_err;
#define mpi_root 0
// end globals

template <class GEN>
vec3<float> nrand_hemisphere(GEN &gen_rad, GEN &gen_phi) {
    float rad = sqrtf(gen_rad());
    float phi = gen_phi() * 2 * M_PI;
    return vec3<float>(rad * cosf(phi), rad * sinf(phi), sqrtf(1 - rad * rad));
}

template <class GEN>
vec3<double> nrand_hemisphere_spherical(GEN &gen_rad, GEN &gen_phi) {
    double theta = acos(sqrt(gen_rad()));
    double phi = gen_phi() * 2 * M_PI;
    return vec3<double>(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}

void init_it(int *argc, char ***argv) {
    mpi_err = MPI_Init(argc, argv);
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

/**
 * Parallel plates
 * param a - x dimension
 * param b - y dimension
 * param c - distance between plates
 */
double real_result(double a, double b, double c) {
    double x = a/c, y = b/c, xq = x * x, yq = y * y;
    double result = 2 / M_PI / x / y * (log(sqrt((1 + xq) * (1 + yq) / (1 + xq + yq))) + x * sqrt(1 + yq) * atan(x / sqrt(1 + yq)) + y * sqrt(1 + xq) * atan(y / sqrt(1 + xq)) - x * atan(x) - y * atan(y));
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
    const int NUM_RAYS = 100, STEPS = 100;
    float c = 1, a = 1, b = 1, a_step = a / STEPS, b_step = b / STEPS;

    auto A1 = poly<float, 4>(0, 0, 0, a, 0, 0, a, b, 0, 0, b, 0);
    auto A2 = poly<float, 4>(0, 0, c, 0, b, c, a, b, c, a, 0, c);
    auto ray = ray3<float>(a / 2, b / 2, 0, 1, 1, 1);

    vec3<float> intr;
    float dist = A2.intersect_ray(ray, intr);
    printf("%f\n", dist);

    /*
    plane3<float> A1(0, 0, 0, 0, 0, 1);
    plane3<float> A2(0, 0, c / 2, 1, 0, 0);

    std::default_random_engine eng1, eng2(1000);
    std::uniform_real_distribution<float> x_distr(0, 1);

    auto rgen = std::bind(x_distr, eng1);
    auto pgen = std::bind(x_distr, eng2);
    unsigned long long ok = 0, total = 0;

    struct timespec t2, t3;
    clock_gettime(CLOCK_MONOTONIC,  &t2);

    for (auto i = a_step / 2; i < a; i += a_step) {
        for (auto j = b_step / 2; j < b; j += b_step) {
            vec3<float> p0(i, j, 0);
            for (int r = 0; r < NUM_RAYS; ++r) {
                // Generate random ray
                auto u = nrand_hemisphere(rgen, pgen);
                ray3<float> ray(p0, u);
                
                // Calculate intersection with A2 plane
                vec3<float> intr;
                float dist = ray.intersect_plane(A2, intr);

                if (dist > 0 &&
                    intr.z >= 0 && intr.z < c &&
                    intr.y >= 0 && intr.y < b) {
                    ++ok;
                }
                ++total;
            }
        }
    }
    auto result = (float)ok/(float)(total);

    clock_gettime(CLOCK_MONOTONIC,  &t3);
    auto dt1 = (t3.tv_sec - t2.tv_sec) + (float) (t3.tv_nsec - t2.tv_nsec) * 1e-9;

    auto real = real_result_90(a, b, c);
    printf("Result: %f, real: %f, rel: %f, rays: %llu\n", result, real, result / real, total);
    printf("%.3f ms\n", dt1 * 1000);
    */

    return 0;
}
