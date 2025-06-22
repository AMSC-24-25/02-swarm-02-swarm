#include <benchmark/benchmark.h>

#include <random>
#include <chrono>

#include "Algorithm.hpp"
#include "Rosenbrock.hpp"

static void TunnellingOpenMP_Rosenbrock(benchmark::State& state) {
	const size_t dimensions = 2;
	const size_t num_positions = state.range(0);
	const size_t num_threads = state.range(1);
	const size_t max_iterations = 1000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	double sigma_max = 7.0;
	double sigma_min = 1.e-8;
	const double gamma = 0.000001;
	const double beta_adjust_factor = 0.7;
	double beta = 50.0;
	const size_t tunnelling = 6;
	const double beta_tresholding = 0.2;
    const size_t time_step_updating = 2000;
	const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

	for (auto _ : state) {
		const auto start = std::chrono::high_resolution_clock::now();

		const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(dimensions, max_iterations, seed, lower_bound, upper_bound, sigma_max, sigma_min, r, gamma, beta_adjust_factor, false, beta, tunnelling, beta_tresholding, num_positions, time_step_updating, num_threads);
        benchmark::ClobberMemory();

		const auto end = std::chrono::high_resolution_clock::now();
		const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		state.SetIterationTime(elapsed_seconds.count());
		state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
	}
}

BENCHMARK(TunnellingOpenMP_Rosenbrock)->ArgNames({"positions","threads"})->ArgsProduct({benchmark::CreateRange(4, 1024, 2), benchmark::CreateRange(1, 8, 2)})->UseManualTime();
