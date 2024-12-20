#include <vector>
#include <cassert>

#include "Sphere.hpp"

double Sphere::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);
	double sum = 0.0;
	for (double x : position) {
		sum += x * x;
	}
	return sum;
}
