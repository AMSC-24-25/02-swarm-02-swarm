#include <vector>
#include <cassert>
#include <functional>
#include <numeric>

#include "Sphere.hpp"

double Sphere::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);

	return std::transform_reduce(position.begin(), position.end(), 0.0, std::plus<double>{},
								 [](const double x) { return x * x; });
}
