#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "funcs.hpp"
#include "time.hpp"

#include <random>
#include "integrate.hpp"
#include "linalg.hpp"

using namespace integrate;
using namespace linalg;

//  globals
int numnodes, myid, mpi_err;
#define mpi_root 0
// end globals

template <class GEN>
vec3<float> rand_hemisphere(GEN &gen_rad, GEN &gen_phi, float rad) {
    float r = gen_rad();
    float phi = gen_phi();
    return vec3<float>(r * cos(phi), r * sin(phi), sqrt(rad * rad - r * r));
}

void init_it(int *argc, char ***argv) {
    mpi_err = MPI_Init(argc, argv);
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

int main(int argc, char **argv) {
    /*
    std::default_random_engine eng;
    std::uniform_real_distribution<float> phi_distr(0, 2 * M_PI);
    std::uniform_real_distribution<float> radius_distr(0, 1);
    auto rgen = std::bind(radius_distr, eng);
    auto pgen = std::bind(phi_distr, eng);

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

    ray3<float> ray(1, 0, 2, 1, 0, -1);
    plane3<float> plane(1, 0, 0.5, 0, 0, 0.5);

    vec3<float> intr;
    auto dist = ray.intersect_plane(plane, intr);
    if (dist > 0) {
        printf("(%f, %f, %f)\n", intr.x, intr.y, intr.z);
    } else {
        printf("No intersection\n");
    }

    vec3<float> intr2;
    auto sim_dist = plane.intersect_ray(ray, intr2);
    assert(dist == sim_dist);
}
