#include <vector>
#include <random>

#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "DifferentialEvolution.hpp"
#include "Candidate.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(GeneticConvergenceOpenMP, Sphere) {
	const size_t dimensions = 2;
	const size_t num_candidates = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double F = 0.5; //Tipico valore: 0.4 - 1.0
	const double CR = 0.8; //Tipico valore: 0.7 - 0.9

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Sphere>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_differential_evolution(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,s,1,true);
	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}

TEST(GeneticConvergenceOpenMP, EuclideanDistance) {
	const size_t dimensions = 2;
	const size_t num_candidates = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double F = 0.5; //Tipico valore: 0.4 - 1.0
	const double CR = 0.8; //Tipico valore: 0.7 - 0.9

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<EuclideanDistance>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_differential_evolution(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,s,1,true);
	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}

TEST(GeneticConvergenceOpenMP, Rosenbrock) {
	const size_t dimensions = 2;
	const size_t num_candidates = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double F = 0.5; //Tipico valore: 0.4 - 1.0
	const double CR = 0.8; //Tipico valore: 0.7 - 0.9

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Rosenbrock>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_differential_evolution(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,s,1,true);
	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}

TEST(GeneticConvergenceOpenMP, Rastrigin) {
	const size_t dimensions = 2;
	const size_t num_candidates = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double F = 0.5; //Tipico valore: 0.4 - 1.0
	const double CR = 0.8; //Tipico valore: 0.7 - 0.9

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Rastrigin>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_differential_evolution(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,s,1,true);
	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}
