#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <utility>
#include <vector>
#include <memory>

#include "ObjectiveFunction.hpp"

namespace algorithm {

/*
 * Runs a particle swarm optimization (PSO) algorithm to optimize the given objective function.
 *
 * Parameters:
 * - dimensions:      The number of dimensions in the solution space.
 * - num_particles:   The number of particles to be used in the swarm.
 * - max_iterations:  The number of iterations to be performed.
 * - seed:            Seed to be used for the random number generation.
 * - lower_bound:     The lower bound of the solution space for all dimensions.
 * - upper_bound:     The upper bound of the solution space for all dimensions.
 * - func:            A pointer to the objective function to be optimized.
 * - n_threads:       The number of OpenMP threads to be used (when possible).
 * - verbose:         If true, enables logging of the best solution at every iteration.
 *
 * Returns a pair of:
 * - An std::vector of doubles representing the best position found in the search space.
 * - A double representing the value of the objective function at the best position.
 */
std::pair<std::vector<double>, double> run_swarm(const size_t dimensions, const size_t num_particles,
												 const size_t max_iterations, const size_t seed,
												 const double lower_bound, const double upper_bound,
												 const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads,
												 const bool verbose);


std::pair<std::vector<double>, double> run_stochastic_tunnelling(const size_t dimensions,
												 const size_t max_iterations, const size_t seed,
												 const double lower_bound, const double upper_bound, const double sigma_max, const double sigma_min,
												 const  ObjectiveFunction& func, const double gamma,
												 const double beta_adjust_factor, const size_t moving_avg_window, const bool verbose, double beta, const size_t tunnelling);

std::pair<std::vector<double>, double> run_differential_evolution(const size_t dimensions,
													const size_t num_candidates, const double lower_bound,
													const double upper_bound, const size_t seed,
													const size_t max_gen, const double F,
													const double CR, const std::unique_ptr<ObjectiveFunction>& func,
													const size_t n_threads, const bool verbose);





/*
 * Runs a parallel genetic algorithm with OpenMP to optimize the given objective function.
 *
 * Parameters:
 * - dimensions:      The number of dimensions in the solution space.
 * - num_creatures:   The number of creatures in the population.
 * - max_iterations:  The number of iterations (or generations) to be performed.
 * - seed:            Seed to be used for the random number generation.
 * - lower_bound:     The lower bound of the solution space for all dimensions.
 * - upper_bound:     The upper bound of the solution space for all dimensions.
 * - mutation_rate:   The probability of mutation for each creature in the population.
 * - survival_rate:   The proportion of creatures that survive each generation (elitism).
 * - func:            A pointer to the objective function to be optimized.
 * - n_threads:       The number of OpenMP threads to be used (when possible).
 * - verbose:         If true, enables logging of the best solution at every iteration.
 *
 * Returns a pair of:
 * - An std::vector of doubles representing the best position found in the search space.
 * - A double representing the value of the objective function at the best position.
 */
std::pair<std::vector<double>, double> run_genetic_openmp(const size_t dimensions, const size_t num_creatures,
														  const size_t max_iterations, const size_t seed,
														  const double lower_bound, const double upper_bound,
														  const double mutation_rate, const double survival_rate,
														  const std::unique_ptr<ObjectiveFunction>& func,
														  const size_t n_threads, const bool verbose);

#if defined(USE_MPI) && USE_MPI == 1

/*
 * Runs a parallel genetic algorithm with MPI to optimize the given objective function.
 *
 * Parameters:
 * - dimensions:      The number of dimensions in the solution space.
 * - num_creatures:   The number of creatures in the population.
 * - max_iterations:  The number of iterations (or generations) to be performed.
 * - seed:            Seed to be used for the random number generation.
 * - lower_bound:     The lower bound of the solution space for all dimensions.
 * - upper_bound:     The upper bound of the solution space for all dimensions.
 * - mutation_rate:   The probability of mutation for each creature in the population.
 * - survival_rate:   The proportion of creatures that survive each generation (elitism).
 * - func:            A pointer to the objective function to be optimized.
 * - verbose:         If true, enables logging of the best solution at every iteration.
 *
 * Returns a pair of:
 * - An std::vector of doubles representing the best position found in the search space.
 * - A double representing the value of the objective function at the best position.
 *
 * Note:
 * The caller must call MPI_Init() before calling this function and MPI_Finalize() after this function returns.
 */
std::pair<std::vector<double>, double> run_genetic_mpi(const size_t dimensions, const size_t num_creatures,
													   const size_t max_iterations, const size_t seed,
													   const double lower_bound, const double upper_bound,
													   const double mutation_rate, const double survival_rate,
													   const std::unique_ptr<ObjectiveFunction>& func,
													   const bool verbose);
#endif	// USE_MPI

}  // namespace algorithm

#endif	// ALGORITHM_HPP