#include <vector>
#include <random>

#include <gtest/gtest.h>
#include <mpi.h>

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

TEST(DeConvergenceMPI, Sphere) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
            algorithm::run_de_mpi(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,s,false);

    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);
        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; i++) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
        }
    }

}

TEST(DeConvergenceMPI, EuclideanDistance) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const size_t dimensions = 2;
    const size_t num_candidates = 100;
    const size_t max_iterations = 1'000;
    const size_t seed = 42;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const double F = 0.5; //Tipico valore: 0.4 - 1.0
    const double CR = 0.8; //Tipico valore: 0.7 - 0.9

    const std::unique_ptr<ObjectiveFunction> e = std::make_unique<EuclideanDistance>();

    const std::pair<std::vector<double>, double> result =
            algorithm::run_de_mpi(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,e,false);
    if(rank==0) {
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);
        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; i++) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
        }
    }

}

TEST(DeConvergenceMPI, Rosenbrock) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const size_t dimensions = 2;
    const size_t num_candidates = 100;
    const size_t max_iterations = 1'000;
    const size_t seed = 42;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const double F = 0.5; //Tipico valore: 0.4 - 1.0
    const double CR = 0.8; //Tipico valore: 0.7 - 0.9

    const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

    const std::pair<std::vector<double>, double> result =
            algorithm::run_de_mpi(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,r,false);

    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; i++) {
            expected_minimum[i] = std::pow(static_cast<Rosenbrock*>(r.get())->a, static_cast<double>(i + 1));
        }
        for (size_t i = 0; i < dimensions; i++) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
        }
    }
}

TEST(DeConvergenceMPI, Rastrigin) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const size_t dimensions = 2;
    const size_t num_candidates = 100;
    const size_t max_iterations = 1'000;
    const size_t seed = 42;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const double F = 0.5; //Tipico valore: 0.4 - 1.0
    const double CR = 0.8; //Tipico valore: 0.7 - 0.9

    const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rastrigin>();

    const std::pair<std::vector<double>, double> result =
            algorithm::run_de_mpi(dimensions, num_candidates ,lower_bound, upper_bound, seed,max_iterations,F,CR,r,false);

    if(rank==0){
        EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

        std::vector<double> expected_minimum(dimensions, 0.0);
        for (size_t i = 0; i < dimensions; i++) {
            EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
        }
    }
}
