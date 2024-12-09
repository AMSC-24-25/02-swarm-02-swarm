#include <cmath>
#include <vector>

#include "EuclideanDistance.hpp"

double EuclideanDistance::operator()(const std::vector<double>& position) const {
	double sum = 0.0;
	for (double x : position) {
		sum += x * x;
	}
	return std::sqrt(sum);
}
