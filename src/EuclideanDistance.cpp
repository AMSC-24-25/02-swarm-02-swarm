#include <cmath>
#include <cassert>
#include <vector>
#include <functional>
#include <numeric>

#include "EuclideanDistance.hpp"

double EuclideanDistance::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);

	return std::sqrt(std::transform_reduce(position.begin(), position.end(), 0.0, std::plus<double>{},
										   [](const double x) { return x * x; }));
}
