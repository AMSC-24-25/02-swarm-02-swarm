#include <benchmark/benchmark.h>

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "Swarm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"

static void Iteration_10000(benchmark::State& state) {
	// Initialize the swarm once before the benchmark loop

	const int dimensions = 6;
	int num_particles = 10000;
	int num_threads = 1;
	const int max_iterations = 100;

	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; ++i) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rosenbrock r;
	static Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, num_threads);

	for (auto _ : state) {
		// Update inertia weight
		swarm.updateInertia(max_iterations, w_min, w_max);

		// Update particle positions and velocities
		swarm.updateParticles();

		// Evaluate fitness function for all particles and find global best
		swarm.findBestFitness();
	}
}

static void Iteration_100(benchmark::State& state) {
	// Initialize the swarm once before the benchmark loop

	const int dimensions = 6;
	int num_particles = 100;
	int num_threads = 1;
	const int max_iterations = 100;

	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; ++i) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rosenbrock r;
	static Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, num_threads);

	for (auto _ : state) {
		// Update inertia weight
		swarm.updateInertia(max_iterations, w_min, w_max);

		// Update particle positions and velocities
		swarm.updateParticles();

		// Evaluate fitness function for all particles and find global best
		swarm.findBestFitness();
	}
}

BENCHMARK(Iteration_10000);
BENCHMARK(Iteration_100);
BENCHMARK_MAIN();