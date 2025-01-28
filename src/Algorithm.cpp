#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <omp.h>

#if defined(USE_MPI) && USE_MPI == 1
#include <mpi.h>

#include "DistributedGeneticAlgorithm.hpp"
#endif	// USE_MPI

#include "Algorithm.hpp"
#include "Particle.hpp"
#include "Swarm.hpp"
#include "GeneticAlgorithm.hpp"
#include "Creature.hpp"

namespace algorithm {

namespace utils {

void print_point(const size_t dimensions, const std::vector<double>& x, const double f) {
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	std::cout << "  f(";
	for (size_t i = 0; i < dimensions; i++) {
		std::cout << std::scientific << x.at(i);
		if (i < dimensions - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << std::scientific << f << std::endl;
	std::cout.copyfmt(oldState);
}

}  // namespace utils

std::pair<std::vector<double>, double> run_swarm(const size_t dimensions, const size_t num_particles,
												 const size_t max_iterations, const size_t seed,
												 const double lower_bound, const double upper_bound,
												 const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads,
												 const bool verbose) {
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const double c1 = 2.0;
	const double c2 = 2.0;
	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Swarm swarm = Swarm(swarmParticles, lower_bound, upper_bound, c1, c2, w, seed, *func, n_threads);

	const double beginning = omp_get_wtime();

	for (size_t i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();

		if (verbose) {
			std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
			std::cout << "  Current minimum: " << std::endl;
			utils::print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
			std::cout << std::endl;
		}
	}

	const double end = omp_get_wtime();

	if (verbose) {
		std::cout << std::endl;
		std::cout << "Minimum found:" << std::endl;
		utils::print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
		std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
				  << std::endl;
		std::cout << std::endl;
	}

	return {swarm.bestGlobalPosition, swarm.minimum};
}

std::pair<std::vector<double>, double> run_genetic_openmp(const size_t dimensions, const size_t num_creatures,
														  const size_t max_iterations, const size_t seed,
														  const double lower_bound, const double upper_bound,
														  const double mutation_rate, const double survival_rate,
														  const std::unique_ptr<ObjectiveFunction>& func,
														  const size_t n_threads, const bool verbose) {
	std::vector<Creature> creatures;

	{
		std::mt19937 rnd{seed};
		std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
		for (size_t i{0}; i < num_creatures; i++) {
			std::vector<double> tmp(dimensions);
			std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
			creatures.push_back(Creature(tmp));
		}
	}

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, *func, n_threads);

	const double beginning = omp_get_wtime();

	// Initial evaluation
	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();

		if (verbose) {
			std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
			std::cout << "  Current minimum: " << std::endl;
			utils::print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
			std::cout << std::endl;
		}
	}

	const double end = omp_get_wtime();

	if (verbose) {
		std::cout << std::endl;
		std::cout << "Minimum found:" << std::endl;
		utils::print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
		std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
				  << std::endl;
		std::cout << std::endl;
	}

	return {ga.bestCreature.position, ga.bestCreature.fitness};
}

#if defined(USE_MPI) && USE_MPI == 1
std::pair<std::vector<double>, double> run_genetic_mpi(const size_t dimensions, const size_t num_creatures,
													   const size_t max_iterations, const size_t seed,
													   const double lower_bound, const double upper_bound,
													   const double mutation_rate, const double survival_rate,
													   const std::unique_ptr<ObjectiveFunction>& func,
													   const bool verbose) {
	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	const size_t local_size = num_creatures / static_cast<size_t>(world_size);

	std::vector<std::vector<double>> creature_positions;

	// Only the "root process" initializes the creatures
	if (world_rank == 0) {
		std::mt19937 rnd{seed};
		std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
		for (size_t i{0}; i < num_creatures; i++) {
			std::vector<double> tmp(dimensions);
			std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
			creature_positions.push_back(tmp);
		}
	}

	// Split the creatures data among processes
	if (world_rank == 0) {
		// The root process sends to everyone
		for (int i = 1; i < world_size; i++) {
			for (size_t j{i * local_size}; j < (i + 1) * local_size; j++) {
				MPI_Send(creature_positions[j].data(), dimensions, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
	} else {
		// The other processes only receive
		creature_positions.resize(local_size);
		for (size_t i{0}; i < local_size; i++) {
			creature_positions[i].resize(dimensions);
			MPI_Recv(creature_positions[i].data(), dimensions, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	// The root process resizes its local vector
	if (world_rank == 0) {
		creature_positions.resize(local_size);
	}

	DistributedGeneticAlgorithm ga(world_rank, world_size, num_creatures, creature_positions, lower_bound, upper_bound,
								   mutation_rate, survival_rate, *func);

	MPI_Barrier(MPI_COMM_WORLD);
	const double beginning = MPI_Wtime();

	// Initial evaluation
	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();

		if (verbose && world_rank == 0) {
			std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
			std::cout << "  Current minimum: " << std::endl;
			utils::print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
							   ga.creature_fitnesses.at(ga.best_creature_index));
			std::cout << std::endl;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	const double end = MPI_Wtime();

	MPI_Finalize();

	if (verbose && world_rank == 0) {
		std::cout << std::endl;
		std::cout << "Minimum found:" << std::endl;
		utils::print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
						   ga.creature_fitnesses.at(ga.best_creature_index));
		std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
				  << std::endl;
		std::cout << std::endl;
	}

	return {ga.creature_positions.at(ga.best_creature_index), ga.creature_fitnesses.at(ga.best_creature_index)};
}
#endif	// USE_MPI

}  // namespace algorithm