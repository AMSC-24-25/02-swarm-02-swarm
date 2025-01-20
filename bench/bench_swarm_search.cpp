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

static void SwarmSearch_Rosenbrock(benchmark::State& state) {
	const int dimensions = 6;
	const int num_particles = state.range(0);
	const int max_iterations = 100;

	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	std::generate_n(std::back_inserter(swarmParticles), num_particles,
					[&]() { return Particle(dimensions, lowerBound, upperBound, 42); });

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rosenbrock r;
	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, 1);

	for (auto _ : state) {
		const auto start = std::chrono::high_resolution_clock::now();

		// Update inertia weight
		swarm.updateInertia(max_iterations, w_min, w_max);
		// Update particle positions and velocities
		swarm.updateParticles();
		// Evaluate fitness function for all particles and find global best
		swarm.findBestFitness();

		const auto end = std::chrono::high_resolution_clock::now();
		const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		state.SetIterationTime(elapsed_seconds.count());
	}
}

BENCHMARK(SwarmSearch_Rosenbrock)->Range(1, 100'000)->RangeMultiplier(10)->UseManualTime();

BENCHMARK_MAIN();