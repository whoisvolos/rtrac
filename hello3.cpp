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

int main(int argc, char **argv) {
    /*

    struct timespec t2, t3;
    clock_gettime(CLOCK_MONOTONIC,  &t2);
    vec3<float> acc(0, 0, 0);
    for (int i = 0; i < 1000000; i++) {
        acc = acc + rand_hemisphere(rgen, pgen, 1);
    }
    clock_gettime(CLOCK_MONOTONIC,  &t3);
    double dt1 = (t3.tv_sec - t2.tv_sec) + (float) (t3.tv_nsec - t2.tv_nsec) * 1e-9;
    printf("%.3f ms\n", dt1 * 1000);
    printf("(%f, %f, %f)\n", acc.x, acc.y, acc.z);
    */

    const int NUM_RAYS = 1000;
    float c = 1, a = 1, b = 1, a_step = 0.1, b_step = 0.1;
    plane3<float> A1(0, 0, 0, 0, 0, 1);
    plane3<float> A2(0, 0, c, 0, 0, -1);

    std::default_random_engine eng;
    std::uniform_real_distribution<float> phi_distr(0, 2 * M_PI);
    std::uniform_real_distribution<float> radius_distr(0, 1);
    auto rgen = std::bind(radius_distr, eng);
    auto pgen = std::bind(phi_distr, eng);

    for (float i = 0; i < a; i += a_step) {
        for (float j = 0; j < b; j += b_step) {
            vec3<float> p0(i + a_step / 2, j + b_step / 2, 0);
            for (int r = 0; r < NUM_RAYS; ++r) {
                // Generate random ray
                vec3<float> u = nrand_hemisphere(rgen, pgen);
                ray3<float> ray(p0, u);
                
                // Calculate intersection with A2 plane
                vec3<float> intr;
                float dist = ray.intersect_plane(A2, intr);
                
            }
//          printf("(%f, %f, %f) (%f, %f, %f), %f\n", p0.x, p0.y, p0.z, u.x, u.y, u.z, u.norm);
        }
    }
}
