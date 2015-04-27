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
    float r = gen_rad();
    float phi = gen_phi();
    return vec3<float>(r * cos(phi), r * sin(phi), sqrt(1 - r * r));
}

void init_it(int *argc, char ***argv) {
    mpi_err = MPI_Init(argc, argv);
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

float real_result(float a, float b, float c) {
    float x = a/c, y = b/c, xq = x * x, yq = y * y;
    float result = 2 / M_PI / x / y * (log(sqrt((1 + xq) * (1 + yq) / (1 + xq + yq))) + x * sqrt(1 + yq) * atan(x / sqrt(1 + yq)) + y * sqrt(1 + xq) * atan(y / sqrt(1 + xq)) - x * atan(x) - y * atan(y));
    return result;
}

int main(int argc, char **argv) {
    const int NUM_RAYS = 100;
    float c = 1, a = 1, b = 1, a_step = 0.01, b_step = 0.01;
    plane3<float> A1(0, 0, 0, 0, 0, 1);
    plane3<float> A2(0, 0, c, 0, 0, -1);

    std::default_random_engine eng;
    std::uniform_real_distribution<float> phi_distr(0, 2 * M_PI);
    std::uniform_real_distribution<float> radius_distr(0, 1);

    std::normal_distribution<float> n_distr(0, 1);

    auto rgen = std::bind(radius_distr, eng);
    auto pgen = std::bind(phi_distr, eng);
    unsigned long long ok = 0, total = 0;

    struct timespec t2, t3;
    clock_gettime(CLOCK_MONOTONIC,  &t2);

    for (float i = a_step / 2; i < a; i += a_step) {
        for (float j = b_step / 2; j < b; j += b_step) {
            vec3<float> p0(i, j, 0);
            for (int r = 0; r < NUM_RAYS; ++r) {
                // Generate random ray
                vec3<float> u = nrand_hemisphere(rgen, pgen);
                ray3<float> ray(p0, u);
                
                // Calculate intersection with A2 plane
                vec3<float> intr;
                float dist = ray.intersect_plane(A2, intr);
                if (dist > 0 &&
                    intr.x >= 0 && intr.x < a &&
                    intr.y >= 0 && intr.y < b) {
                    ++ok;
                }
                ++total;
            }
        }
    }
    float result = (float)ok/(float)(total);

    clock_gettime(CLOCK_MONOTONIC,  &t3);
    double dt1 = (t3.tv_sec - t2.tv_sec) + (float) (t3.tv_nsec - t2.tv_nsec) * 1e-9;

    float real = real_result(a, b, c);
    printf("Result: %f, real: %f, rel: %f, rays: %llu\n", result, real, result / real, total);
    printf("%.3f ms\n", dt1 * 1000);
}
