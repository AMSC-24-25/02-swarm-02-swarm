#include <vector>
#include <random>
#include <memory>
#include <cmath>
#include <cassert>
#include <mpi.h>
#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "SimulatedAnnealing.hpp"
#include "State.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
    assert(std::isfinite(expected));
    assert(std::isfinite(actual));
    return std::abs(expected - actual);
}

TEST(SaConvergenceMPI, Sphere) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const size_t dimensions = 2;
    const size_t max_iterations = 1000;
    const size_t dwell = 200;
    const double initial_temperature = 8.0;
    const double temperature_scale = 0.95;
    const double initial_step_size = 2.0;
    const double step_size_scale = 0.98;
    const double boltzmann_k = 1.0;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const size_t seed = 42;

    const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Sphere>();

    const std::pair<std::vector<double>, double> result =
        algorithm::run_sa_mpi(
            dimensions, max_iterations, dwell,
            initial_temperature, temperature_scale,
            initial_step_size, step_size_scale,
            boltzmann_k, lower_bound, upper_bound,
            s, seed, false
        );
    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);
        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; ++i) {
        EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
        }
    }    
}

TEST(SaConvergenceMPI, EuclideanDistance) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const size_t dimensions = 2;
    const size_t max_iterations = 3000;
    const size_t dwell = 500;
    const double initial_temperature = 10.0;
    const double temperature_scale = 0.99;
    const double initial_step_size = 0.1;
    const double step_size_scale = 0.99;
    const double boltzmann_k = 1.0;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const size_t seed = 42;

    const std::unique_ptr<ObjectiveFunction> e = std::make_unique<EuclideanDistance>();

    const std::pair<std::vector<double>, double> result =
        algorithm::run_sa_mpi(
            dimensions, max_iterations, dwell,
            initial_temperature, temperature_scale,
            initial_step_size, step_size_scale,
            boltzmann_k, lower_bound, upper_bound,
            e, seed, false
        );

    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 2e-3);

        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; ++i) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 2e-3);
        }
    }
}

TEST(SaConvergenceMPI, Rosenbrock) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const size_t dimensions = 2;
    const size_t max_iterations = 1000;
    const size_t dwell = 200;
    const double initial_temperature = 15.0;
    const double temperature_scale = 0.93;
    const double initial_step_size = 0.5;
    const double step_size_scale = 0.99;
    const double boltzmann_k = 1.0;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const size_t seed = 42;
   

    const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

    const std::pair<std::vector<double>, double> result =
        algorithm::run_sa_mpi(
            dimensions, max_iterations, dwell,
            initial_temperature, temperature_scale,
            initial_step_size, step_size_scale,
            boltzmann_k, lower_bound, upper_bound,
            r, seed, false
        );

    if(rank==0){    
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; ++i) {
            expected_minimum[i] = std::pow(static_cast<Rosenbrock*>(r.get())->a, i + 1);
        }
        for (size_t i = 0; i < dimensions; ++i) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-2);
        }
    }
}

TEST(SaConvergenceMPI, Rastrigin) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const size_t dimensions = 2;
    const size_t max_iterations = 5000;
    const size_t dwell = 500;
    const double initial_temperature = 20.0;
    const double temperature_scale = 0.97;
    const double initial_step_size = 0.3;
    const double step_size_scale = 0.99;
    const double boltzmann_k = 1.0;
    const double lower_bound = -5.12;
    const double upper_bound = 5.12;
    const size_t seed = 42;


    const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rastrigin>();

    const std::pair<std::vector<double>, double> result =
        algorithm::run_sa_mpi(
            dimensions, max_iterations, dwell,
            initial_temperature, temperature_scale,
            initial_step_size, step_size_scale,
            boltzmann_k, lower_bound, upper_bound,
            r, seed, false
        );
    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 3e-3);

        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; ++i) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 3e-3);
        }
    }
}

