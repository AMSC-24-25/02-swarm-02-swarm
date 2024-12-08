#include <vector>

#include "ObjectiveFunction.hpp"

double ObjectiveFunction::operator()(const std::vector<double>& position) {
	double sum = 0.0;
	for (double x : position) {
		sum += x * x;
	}
	return sum;
}
