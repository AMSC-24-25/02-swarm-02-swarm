#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

#include "Swarm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(FunctionConvergence, Sphere) {
	const int dimensions = 2;
	const int num_particles = 100;
	const int max_iterations = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Sphere s;
	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, s, 1);

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
	}

	EXPECT_LE(absolute_error(0.0, swarm.minimum), 1e-9);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (int i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), swarm.bestGlobalPosition.at(i)), 1e-9);
	}
}

TEST(FunctionConvergence, EuclideanDistance) {
	const int dimensions = 2;
	const int num_particles = 100;
	const int max_iterations = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	EuclideanDistance ed;
	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, ed, 1);

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
	}

	EXPECT_LE(absolute_error(0.0, swarm.minimum), 1e-9);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (int i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), swarm.bestGlobalPosition.at(i)), 1e-9);
	}
}

TEST(FunctionConvergence, Rosenbrock) {
	const int dimensions = 2;
	const int num_particles = 100;
	const int max_iterations = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rosenbrock r;
	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, 1);

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
	}

	EXPECT_LE(absolute_error(0.0, swarm.minimum), 1e-6);

	std::vector<double> expected_minimum(dimensions);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = std::pow(r.a, static_cast<double>(i + 1));
	}
	for (int i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), swarm.bestGlobalPosition.at(i)), 1e-3);
	}
}

TEST(FunctionConvergence, Rastrigin) {
	const int dimensions = 2;
	const int num_particles = 100;
	const int max_iterations = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rastrigin r;
	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, 1);

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
	}

	EXPECT_LE(absolute_error(0.0, swarm.minimum), 1e-8);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (int i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), swarm.bestGlobalPosition.at(i)), 1e-8);
	}
}
