#include <vector>
#include <cmath>

#include <mpi.h>

#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "Sphere.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(GeneticConvergenceMPI, Sphere) {

	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t max_iterations = 1000;
	double sigma_max = 1.0;
	double sigma_min = 1.e-5;
	const double gamma = 0.00005;
	const double beta_adjust_factor = 0.8;
	double beta = 10000.0;
	const size_t tunnelling = 12;
	const double beta_tresholding = 0.1;
    const size_t num_positions = 100;
    const size_t time_step_updating = 200;

	const ObjectiveFunction& s = Sphere();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_tunnelling_mpi(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, true,beta, tunnelling, beta_tresholding, num_positions, time_step_updating);


	EXPECT_LE(absolute_error(0.0, result.second), 1e-1);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
	}

}


TEST(GeneticConvergenceMPI, Rosenbrock) {
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 36;
	const size_t max_iterations = 1000;
	double sigma_max = 1.0;
	double sigma_min = 1.e-4;
	const double gamma = 0.0001;
	const double beta_adjust_factor = 0.9;
	double beta = 500.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.1;
    const size_t num_positions = 100;
    const size_t time_step_updating = 150;

	const ObjectiveFunction& s = Rosenbrock();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_tunnelling_mpi(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, true,beta, tunnelling, beta_tresholding, num_positions, time_step_updating);


	EXPECT_LE(absolute_error(0.0, result.second), 1e-2);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = 1;
	}
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
	}


}



TEST(GeneticConvergenceMPI, Rastrigin) {
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const size_t dimensions = 2;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t seed = 42;
	const size_t max_iterations = 1000;
	double sigma_max = 1.0;
	double sigma_min = 1.e-4;
	const double gamma = 0.0001;
	const double beta_adjust_factor = 0.9;
	double beta = 500.0;
	const size_t tunnelling = 10;
	const double beta_tresholding = 0.2;
    const size_t num_positions = 100;
    const size_t time_step_updating = 100;

	const ObjectiveFunction& s = Rastrigin();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_tunnelling_mpi(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, s, gamma, beta_adjust_factor, true,beta, tunnelling, beta_tresholding, num_positions, time_step_updating);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-1);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
	}
}