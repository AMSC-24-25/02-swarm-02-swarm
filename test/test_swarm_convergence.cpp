#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

#include "Algorithm.hpp"
#include "Particle.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(SwarmConvergence, Sphere) {
	const size_t dimensions = 2;
	const size_t num_particles = 100;
	const size_t max_iterations = 200;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const std::unique_ptr<ObjectiveFunction> s = std::make_unique<Sphere>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_swarm(dimensions, num_particles, max_iterations, seed, lower_bound, upper_bound, s, 1);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}

TEST(SwarmConvergence, EuclideanDistance) {
	const size_t dimensions = 2;
	const size_t num_particles = 100;
	const size_t max_iterations = 200;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const std::unique_ptr<ObjectiveFunction> ed = std::make_unique<EuclideanDistance>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_swarm(dimensions, num_particles, max_iterations, seed, lower_bound, upper_bound, ed, 1);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}

TEST(SwarmConvergence, Rosenbrock) {
	const size_t dimensions = 2;
	const size_t num_particles = 100;
	const size_t max_iterations = 200;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_swarm(dimensions, num_particles, max_iterations, seed, lower_bound, upper_bound, r, 1);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = std::pow(static_cast<Rosenbrock*>(r.get())->a, static_cast<double>(i + 1));
	}
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 0.1);
	}
}

TEST(SwarmConvergence, Rastrigin) {
	const size_t dimensions = 2;
	const size_t num_particles = 100;
	const size_t max_iterations = 200;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const std::unique_ptr<ObjectiveFunction> ed = std::make_unique<EuclideanDistance>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_swarm(dimensions, num_particles, max_iterations, seed, lower_bound, upper_bound, ed, 1);

	EXPECT_LE(absolute_error(0.0, result.second), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), result.first.at(i)), 1e-3);
	}
}
