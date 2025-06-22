#include <vector>
#include <cmath>

#include <mpi.h>

#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(GeneticConvergenceMPI, Sphere) {
	MPI_Init(NULL, NULL);

	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Sphere>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_genetic_mpi(dimensions, num_creatures, max_iterations, seed, lower_bound, upper_bound,
								   mutation_rate, survival_rate, s, false);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}

	MPI_Finalize();
}

TEST(GeneticConvergenceMPI, EuclideanDistance) {
	MPI_Init(NULL, NULL);

	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	const std::unique_ptr<ObjectiveFunction> ed = std::make_unique<EuclideanDistance>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_genetic_mpi(dimensions, num_creatures, max_iterations, seed, lower_bound, upper_bound,
								   mutation_rate, survival_rate, ed, false);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}

	MPI_Finalize();
}

TEST(GeneticConvergenceMPI, Rosenbrock) {
	MPI_Init(NULL, NULL);

	const size_t dimensions = 2;
	const size_t num_creatures = 1'000;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_genetic_mpi(dimensions, num_creatures, max_iterations, seed, lower_bound, upper_bound,
								   mutation_rate, survival_rate, r, false);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = std::pow(static_cast<Rosenbrock*>(r.get())->a, static_cast<double>(i + 1));
	}
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
	}

	MPI_Finalize();
}

TEST(GeneticConvergenceMPI, Rastrigin) {
	MPI_Init(NULL, NULL);

	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rastrigin>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_genetic_mpi(dimensions, num_creatures, max_iterations, seed, lower_bound, upper_bound,
								   mutation_rate, survival_rate, r, false);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}

	MPI_Finalize();
}
