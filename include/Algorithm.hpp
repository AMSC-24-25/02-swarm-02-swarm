#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <memory>

#include "ObjectiveFunction.hpp"

namespace algorithm {

void run_swarm(const size_t dimensions, const size_t num_particles, const size_t max_iterations, const size_t seed,
			   const double lower_bound, const double upper_bound, const std::unique_ptr<ObjectiveFunction>& func,
			   const size_t n_threads);

void run_genetic_openmp(const size_t dimensions, const size_t num_creatures, const size_t max_iterations,
						const size_t seed, const double lower_bound, const double upper_bound,
						const double mutation_rate, const double survival_rate,
						const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads);

#if defined(USE_MPI) && USE_MPI == 1
std::pair<std::vector<double>, double> run_genetic_mpi(const size_t dimensions, const size_t num_creatures,
													   const size_t max_iterations, const size_t seed,
													   const double lower_bound, const double upper_bound,
													   const double mutation_rate, const double survival_rate,
													   const std::unique_ptr<ObjectiveFunction>& func);
#endif	// USE_MPI

}  // namespace algorithm

#endif	// ALGORITHM_HPP