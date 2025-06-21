
#include <benchmark/benchmark.h>
#include <random>
#include <chrono>
#include "Algorithm.hpp"
#include "Rosenbrock.hpp"  // Cambia funzione qui se vuoi

static void FireflyBFGS_Rosenbrock(benchmark::State& state) {
	const size_t dimensions = 2;
	const size_t num_fireflies = state.range(0);    // Primo range: fireflies
	const size_t num_threads  = state.range(1);     // Secondo range: threads
	const size_t max_iterations = 100;
	const size_t seed = 42;
	const double lower_bound = -5.0;
	const double upper_bound = 5.0;

	// Parametri firefly
	const double alpha = 0.4;
	const double beta  = 7.0;
	const double gamma = 1.0;
	const bool use_cuda = false;

	const std::unique_ptr<ObjectiveFunction> func = std::make_unique<Rosenbrock>();

	for (auto _ : state) {
		const auto start = std::chrono::high_resolution_clock::now();

		auto result = algorithm::run_firefly_bfgs(
			dimensions, num_fireflies, max_iterations,
			seed, lower_bound, upper_bound,
			func, num_threads, false, // verbose = false
			use_cuda, alpha, beta, gamma
		);

		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();

		const auto end = std::chrono::high_resolution_clock::now();
		const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		state.SetIterationTime(elapsed_seconds.count());
		state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
	}
}

BENCHMARK(FireflyBFGS_Rosenbrock)
	->ArgNames({"fireflies","threads"})
	->Ranges({{4, 200},    // fireflies: da 4 a 200, raddoppia
			  {1, 8}})     // threads: 1,2,4,8,16
	->RangeMultiplier(2)
	->UseManualTime();
