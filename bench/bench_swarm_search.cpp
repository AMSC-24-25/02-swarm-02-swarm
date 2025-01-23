#include <benchmark/benchmark.h>

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "Swarm.hpp"
#include "Rosenbrock.hpp"

static void SwarmSearch_Rosenbrock(benchmark::State& state) {
	const int dimensions = 6;
	const int num_particles = state.range(0);
	const int max_iterations = 100;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	const size_t seed = 42;

	for (auto _ : state) {
		std::vector<Particle> swarmParticles;

		for (size_t i = 0; i < num_particles; i++) {
			swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
		}

		const double w_max = 0.9;
		const double w_min = 0.4;
		const double w = w_max;

		Rosenbrock r;
		Swarm swarm = Swarm(swarmParticles, lower_bound, upper_bound, 2.0, 2.0, w, seed, r, 1);

		const auto start = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < max_iterations; i++) {
			swarm.updateInertia(max_iterations, w_min, w_max);
			swarm.updateParticles();
			swarm.findBestFitness();
		}
		benchmark::DoNotOptimize(swarm.minimum);
		benchmark::DoNotOptimize(swarm.bestGlobalPosition);
		benchmark::ClobberMemory();

		const auto end = std::chrono::high_resolution_clock::now();
		const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		state.SetIterationTime(elapsed_seconds.count());
		state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
	}
}

BENCHMARK(SwarmSearch_Rosenbrock)->Range(1, 1'000)->RangeMultiplier(2)->UseManualTime();
