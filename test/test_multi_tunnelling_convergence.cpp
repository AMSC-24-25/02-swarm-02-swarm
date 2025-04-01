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
	const size_t seed = 36;
	const size_t max_iterations = 1000;
	double sigma_max = 10.0;
	double sigma_min = 1.e-9;
	const double gamma = 0.0000005;
	const double beta_adjust_factor = 0.8;
	double beta = 1000.0;
	const size_t tunnelling = 12;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 10;

	const ObjectiveFunction& s = Sphere();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, true,beta, tunnelling, beta_tresholding, num_positions, time_step_updating);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
	}


}

TEST(TunnellingConvergence, Rosenbrock){
	const size_t dimensions = 2;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 40;
	const size_t max_iterations = 1000;
	double sigma_max = 12.0;
	double sigma_min = 1.e-9;
	const double gamma = 0.0000005;
	const double beta_adjust_factor = 0.9;
	double beta = 10000.0;
	const size_t tunnelling = 12;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 10;

	const ObjectiveFunction& s = Rosenbrock();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, true,beta, tunnelling, beta_tresholding, num_positions, time_step_updating);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
	}


}