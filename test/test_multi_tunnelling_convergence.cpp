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
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t max_iterations = 1000;
	const size_t num_threads = 4;
	double sigma_max = 1.0;
	double sigma_min = 1.e-5;
	const double gamma = 0.00005;
	const double beta_adjust_factor = 0.8;
	double beta = 10000.0;
	const size_t tunnelling = 12;
	const double beta_tresholding = 0.1;
    const size_t num_positions = 100;
    const size_t time_step_updating = 200;

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Sphere>();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, false, beta, tunnelling, beta_tresholding, num_positions, time_step_updating, num_threads);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);

	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.5);
	}



}

TEST(TunnellingConvergence, Rosenbrock){
	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t max_iterations = 1000;
	const size_t num_threads = 4;
	double sigma_max = 1.0;
	double sigma_min = 1.e-4;
	const double gamma = 0.0001;
	const double beta_adjust_factor = 0.9;
	double beta = 500.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 150;

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Rosenbrock>();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, false, beta, tunnelling, beta_tresholding, num_positions, time_step_updating, num_threads);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = std::pow(static_cast<Rosenbrock*>(s.get())->a, static_cast<double>(i + 1));
	}
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.5);
	}



}


TEST(TunnellingConvergence, Rastrigin){
	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t num_threads = 4;
	const size_t max_iterations = 1000;
	double sigma_max = 1.0;
	double sigma_min = 5.e-5;
	const double gamma = 0.0001;
	const double beta_adjust_factor = 0.9;
	double beta = 500.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 100;

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Rastrigin>();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, false, beta, tunnelling, beta_tresholding, num_positions, time_step_updating, num_threads);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-1);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.5);
	}



}



TEST(TunnellingConvergence, EuclideanDistance){
	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t num_threads = 4;
	const size_t max_iterations = 1000;
	double sigma_max = 1.0;
	double sigma_min = 5.e-5;
	const double gamma = 0.0001;
	const double beta_adjust_factor = 0.9;
	double beta = 500.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 100;

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<EuclideanDistance>();


	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, false, beta, tunnelling, beta_tresholding, num_positions, time_step_updating, num_threads);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.5);
	}



}