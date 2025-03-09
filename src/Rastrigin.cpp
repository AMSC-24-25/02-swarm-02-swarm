#include <cmath>
#include <cassert>
#include <vector>
#include <functional>
#include <numeric>

#include "Rastrigin.hpp"

/*
 * Multi-dimensional generalization of Rastrigin function.
 *
 * Global minimum at f(0, 0, ..., 0) = 0.0.
 *
 */
double Rastrigin::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);

	return std::transform_reduce(position.begin(), position.end(), A * position.size(), std::plus<double>{},
								 [&A = this->A](const double x) { return x * x - A * std::cos(2.0 * M_PI * x); });
}
