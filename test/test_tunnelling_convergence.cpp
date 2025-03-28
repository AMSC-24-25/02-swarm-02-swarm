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
	const size_t max_iterations = 700;
	double sigma_max = 5.0;
	double sigma_min = 1.e-6;
	const double gamma = 0.000005;
	const double beta_adjust_factor = 0.9;
	const size_t moving_avg_window = 30;
	double beta = 1000.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.3;

	const ObjectiveFunction& s = Sphere();

	const std::pair<std::vector<double>, double> result = 
		algorithm::run_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, moving_avg_window, true,beta, tunnelling, beta_tresholding);

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
	const size_t seed = 28;
	const size_t max_iterations = 1000;
	double sigma_max = 2.0;
	double sigma_min = 1.e-8;
	const double gamma = 0.00001;
	const double beta_adjust_factor = 0.8;
	const size_t moving_avg_window = 30;
	double beta = 5000.0;
	const size_t tunnelling = 12;
	const double beta_thresholding = 0.28;

	const ObjectiveFunction& s = Rosenbrock();

	const std::pair<std::vector<double>, double> result = 
		algorithm::run_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, moving_avg_window, true, beta, tunnelling, beta_thresholding);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
	}


}


/*TEST(TunnellingConvergence, Rastrigin){
	const size_t dimensions = 2;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 36;
	const size_t max_iterations = 1000;
	double sigma_max = 2.0;
	double sigma_min = 1.e-8;
	const double gamma = 0.00001;
	const double beta_adjust_factor = 0.8;
	const size_t moving_avg_window = 15;
	double beta = 5000.0;
	const size_t tunnelling = 10;
	const double beta_thresholding = 0.2;

	const ObjectiveFunction& s = Rastrigin();

	const std::pair<std::vector<double>, double> result = 
		algorithm::run_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, moving_avg_window, true, beta, tunnelling, beta_thresholding);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
	}


}*/