#include <math.h>
#include <random>
#include <assert.h>

namespace integrate {
	typedef double (*integral_func)(double);
	double monte_carlo(integral_func i_func, double x_from, double x_to, double max_y, int iterations);
	double trapezoid(integral_func f, double x_from, double x_to, int iterations);
	double simpson(integral_func f, double x_from, double x_to, int iterations);
}