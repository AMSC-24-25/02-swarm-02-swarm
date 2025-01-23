#include <benchmark/benchmark.h>

#include <random>
#include <chrono>

#include "GeneticAlgorithm.hpp"
#include "Rosenbrock.hpp"

static void Genetic_Rosenbrock(benchmark::State& state) {
	const size_t dimensions = 2;
	const size_t num_creatures = state.range(0);
	const size_t max_iterations = 100;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;

	for (auto _ : state) {
		std::vector<Creature> creatures;

		std::mt19937 rnd{seed};
		std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
		for (size_t i{0}; i < num_creatures; i++) {
			std::vector<double> tmp(dimensions);
			std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
			creatures.push_back(Creature(tmp));
		}

		const double mutation_rate = 0.2;
		const double survival_rate = 0.5;

		Rosenbrock r;

		GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, r, 1);

		const auto start = std::chrono::high_resolution_clock::now();

		ga.evaluateCreatures();
		ga.sortCreatures();

		for (size_t i{0}; i < max_iterations; i++) {
			ga.applyCrossover(seed + i);
			ga.applyMutation(seed + i + 1);
			ga.evaluateCreatures();
			ga.sortCreatures();
		}
		benchmark::DoNotOptimize(ga.bestCreature);
		benchmark::ClobberMemory();

		const auto end = std::chrono::high_resolution_clock::now();
		const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		state.SetIterationTime(elapsed_seconds.count());
		state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
	}
}

BENCHMARK(Genetic_Rosenbrock)->Range(3, 1'000)->RangeMultiplier(2)->UseManualTime();
