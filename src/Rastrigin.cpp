#include <cmath>
#include <cassert>
#include <vector>

#include "Rastrigin.hpp"

/*
 * Multi-dimensional generalization of Rastrigin function.
 *
 * Global minimum at f(0, 0, ..., 0) = 0.0.
 *
 */
double Rastrigin::operator()(const std::vector<double>& position) const {
	assert(position.size() > 0);

	double s = A * position.size();	 // Usa A definita nella classe
	for (const double x : position) {
		s += x * x - A * std::cos(2 * M_PI * x);
	}
	return s;
}
