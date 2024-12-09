#include <vector>

#include "Sphere.hpp"

double Sphere::operator()(const std::vector<double>& position) const {
	double sum = 0.0;
	for (double x : position) {
		sum += x * x;
	}
	return sum;
}
