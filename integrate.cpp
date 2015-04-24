#include "integrate.hpp"

double integrate::monte_carlo(integrate::integral_func i_func, double x_from, double x_to, double max_y, int iterations) {
    assert(x_to > x_from);
    double square = max_y * (x_to - x_from);

    std::default_random_engine eng;
    std::uniform_real_distribution<double> dist_x(x_from, x_to);
    std::uniform_real_distribution<double> dist_y(0.0, max_y);

    double under = 0, above = 0;    

    for (int i = 0; i < iterations; ++i) {
        double x = dist_x(eng), y = dist_y(eng);
        double f_val = i_func(x);
        if (y <= f_val) {
            ++under;
        } else {
            ++above;
        }
    }

    return square * ((double)under / (double)(above + under));
}

double integrate::trapezoid(integrate::integral_func f, double x_from, double x_to, int iterations) {
    assert(x_to > x_from);
    double x_cur = x_from, eps = (x_to - x_from) / (double)iterations;
    double a = f(x_cur);
    a = a != a ? 0.0 : a;
    double result = 0.0;
    for (int i = 0; i < iterations; ++i) {
        x_cur += eps;
        double b = f(x_cur);
        b = b != b ? 0.0 : b;
        result += eps * (a + b) / 2;
        a = b;
    }
    return result;
}

double integrate::simpson(integrate::integral_func f, double x_from, double x_to, int iterations) {
    assert(iterations % 2 == 0);
    assert(x_to > x_from);
    double eps = (x_from - x_to) / iterations;
    double result = 0.0;

    for (int i = 0, n = iterations / 2; i < n; ++i) {
        double a = f(x_from + (i - 1) * 2 * eps);
        double b = f(x_from + (i * 2 - 1) * eps) * 4;
        double c = f(x_from + i * 2 * eps);

        result += (a != a ? 0.0 : a) + (b != b ? 0.0 : b) + (c != c ? 0.0 : c);
    }

    return result * eps / 3;
}