#include <cmath>
#include <cassert>
#include <vector>

#include "EuclideanDistance.hpp"

double EuclideanDistance::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);
	// since you have used other algorithms of the standard library, here you could have used transform_reduce
	double sum = 0.0;
	for (double x : position) {
		sum += x * x;
	}
	return std::sqrt(sum);
}
