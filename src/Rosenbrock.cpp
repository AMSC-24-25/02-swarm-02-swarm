#include <cmath>
#include <cassert>
#include <vector>

#include "Rosenbrock.hpp"

/*
 * Multi-dimensional generalization of Rosenbrock function.
 *
 * Global minimum in f(a; a^2; a^3; a^4; ...) = 0.0.
 *
 * References:
 * https://en.wikipedia.org/wiki/Rosenbrock_function
 * https://www.sfu.ca/~ssurjano/rosen.html
 */
double Rosenbrock::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);

	double s = 0.0;
	for (size_t i = 0; i < position.size() - 1; i++) {
		const double x = position.at(i);
		const double y = position.at(i + 1);
		const double t1 = a - x;
		const double t2 = y - x * x;
		s += t1 * t1 + b * t2 * t2;
	}
	return s;
}
