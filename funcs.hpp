#include <math.h>

inline double planck(double lambda) {
    return pow(lambda, 3) / (exp(lambda) - 1);
}

inline double dome(double x) {
    if (x < -1 || x > 1) {
        return 0.0;
    } else {
        return sqrt(1 - x * x);
    }
}