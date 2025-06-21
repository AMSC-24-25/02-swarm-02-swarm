#include <benchmark/benchmark.h>
#include <chrono>
#include <vector>
#include <memory>

#include "Algorithm.hpp"
#include "Rosenbrock.hpp"

static void SAOpenMP_Rosenbrock(benchmark::State& state) {
    const size_t dimensions = 2;

    const size_t n_threads = state.range(0);
    const size_t max_iterations = state.range(1);

    const size_t dwell_iterations = 300;
    const double initial_temperature = 15.0;
    const double temperature_scale = 0.93;
    const double initial_step_size = 0.5;
    const double step_size_scale = 0.99;
    const double boltzmann_constant = 1.0;

    const std::vector<double> initial_guess(dimensions, -5.0);  
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const size_t seed = 42;
    


    const std::unique_ptr<ObjectiveFunction> func = std::make_unique<Rosenbrock>();

    for (auto _ : state) {
        const auto start = std::chrono::high_resolution_clock::now();

        auto result = algorithm::run_simulated_annealing(
            dimensions,
            max_iterations,
            dwell_iterations,
            initial_temperature,
            temperature_scale,
            initial_step_size,
            step_size_scale,
            boltzmann_constant,
            initial_guess,
            lower_bound,
            upper_bound,
            func,
            seed,
            n_threads,
            false // verbose
        );

        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();

        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
        state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
        
        
        
    }
}

BENCHMARK(SAOpenMP_Rosenbrock)
    ->ArgNames({"threads", "max_iterations"})
    ->Args({1, 200})
    ->Args({2, 200})
    ->Args({4, 200})
    ->Args({8, 200})
    ->Args({1, 512})
    ->Args({2, 512})
    ->Args({4, 512})
    ->Args({8, 512})
    ->Args({1, 1000})
    ->Args({2, 1000})
    ->Args({4, 1000})
    ->Args({8, 1000})
    ->Args({1, 2000})
    ->Args({2, 2000})
    ->Args({4, 2000})
    ->Args({8, 2000})
    ->UseManualTime(); 



