#include <vector>
#include <random>

#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "Position.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(TunnellingConvergence, Sphere){
	const size_t dimensions = 2;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;
	const size_t max_iterations = 200;
	const double f_thresh = 0.3;
	const double sigma = 10.0;
	const double gamma = 0.4;
	const double beta_adjust_factor = 0.9;
	const size_t moving_avg_window = 30;

	const ObjectiveFunction& s = Sphere();

	const std::pair<std::vector<double>, double> result = 
		algorithm::run_stochastic_tunnelling(dimensions, max_iterations, seed, f_thresh, lower_bound, upper_bound, sigma, s, gamma, beta_adjust_factor, moving_avg_window, true);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}


}