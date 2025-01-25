#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <mpi.h>

#include "DistributedGeneticAlgorithm.hpp"

DistributedGeneticAlgorithm::DistributedGeneticAlgorithm(const int world_rank_, const int world_size_,
														 const size_t total_creatures_,
														 const std::vector<std::vector<double>>& creature_positions_,
														 const double lower_bound_, const double upper_bound_,
														 const double survival_rate_, const double mutation_rate_,
														 ObjectiveFunction& func_)
	: world_rank(world_rank_),
	  world_size(world_size_),
	  total_creatures(total_creatures_),
	  creature_positions(creature_positions_),
	  creature_fitnesses(std::vector<double>(creature_positions.size(), std::numeric_limits<double>::infinity())),
	  lower_bound(lower_bound_),
	  upper_bound(upper_bound_),
	  survival_rate(survival_rate_),
	  mutation_rate(mutation_rate_),
	  best_creature_index(0),
	  func(func_) {
	assert(world_rank >= 0);
	assert(world_size >= world_rank);
	assert(total_creatures > 0);
	assert(creature_positions.size() > 0);
	assert(creature_positions.size() < total_creatures);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);
	assert(survival_rate > 0.0 && survival_rate < 1.0);
	assert(mutation_rate >= 0.0 && mutation_rate < 1.0);
}

void DistributedGeneticAlgorithm::evaluateCreatures() {
	// Each process evaluates its own vector of creatures
	for (size_t i{0}; i < creature_positions.size(); i++) {
		if (std::isfinite(creature_fitnesses.at(i))) {
			continue;
		}

		creature_fitnesses[i] = func(creature_positions[i]);
	}
}

void DistributedGeneticAlgorithm::sortCreatures() {
	const size_t dimensions = creature_positions.at(0).size();
	const size_t local_size = total_creatures / static_cast<size_t>(world_size);

	MPI_Barrier(MPI_COMM_WORLD);

	// Each process sends its vector of creatures to the root process
	if (world_rank == 0) {
		creature_positions.resize(total_creatures);
		creature_fitnesses.resize(total_creatures);
		for (int i = 1; i < world_size; i++) {
			for (size_t j{i * local_size}; j < (i + 1) * local_size; j++) {
				creature_positions.at(j).resize(dimensions);
				std::vector<double> tmp = creature_positions[j];
				MPI_Recv(creature_positions.at(j).data(), dimensions, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
						 MPI_STATUS_IGNORE);
			}
			MPI_Recv(creature_fitnesses.data() + (i * local_size), local_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					 MPI_STATUS_IGNORE);
		}
	} else {
		for (size_t i{0}; i < creature_positions.size(); i++) {
			MPI_Send(creature_positions.at(i).data(), dimensions, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Send(creature_fitnesses.data(), creature_positions.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	// The root process then sorts all the creatures
	if (world_rank == 0) {
		std::vector<size_t> indices(total_creatures);
		for (size_t i = 0; i < indices.size(); ++i) {
			indices[i] = i;
		}

		std::sort(indices.begin(), indices.end(),
				  [&](const size_t i, const size_t j) { return creature_fitnesses[i] < creature_fitnesses[j]; });

		std::vector<std::vector<double>> sorted_positions(total_creatures);
		std::vector<double> sorted_fitnesses(total_creatures);
		for (size_t i = 0; i < indices.size(); i++) {
			sorted_positions[i] = creature_positions.at(indices[i]);
			sorted_fitnesses[i] = creature_fitnesses.at(indices[i]);
		}

		creature_positions = sorted_positions;
		creature_fitnesses = sorted_fitnesses;
	}

	// The root process then sends back all the creatures to the processes
	if (world_rank == 0) {
		for (int i = 1; i < world_size; i++) {
			for (size_t j{i * local_size}; j < (i + 1) * local_size; j++) {
				MPI_Send(creature_positions.at(j).data(), dimensions, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(creature_fitnesses.data() + (i * local_size), local_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	} else {
		for (size_t i{0}; i < local_size; i++) {
			creature_positions[i].resize(dimensions);
			MPI_Recv(creature_positions[i].data(), dimensions, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		MPI_Recv(creature_fitnesses.data(), local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (world_rank == 0) {
		creature_positions.resize(local_size);
		creature_fitnesses.resize(local_size);
	}
}

void DistributedGeneticAlgorithm::applyCrossover(const size_t) {
	std::cout << "[" << world_rank << "] TODO: implement crossover" << std::endl;
}

void DistributedGeneticAlgorithm::applyMutation(const size_t) {
	std::cout << "[" << world_rank << "] TODO: implement mutation" << std::endl;
}
