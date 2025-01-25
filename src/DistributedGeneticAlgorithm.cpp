#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>

#include "DistributedGeneticAlgorithm.hpp"

DistributedGeneticAlgorithm::DistributedGeneticAlgorithm(const int world_rank_,
														 const std::vector<std::vector<double>>& creature_positions_,
														 const double lower_bound_, const double upper_bound_,
														 const double survival_rate_, const double mutation_rate_,
														 ObjectiveFunction& func_)
	: world_rank(world_rank_),
	  creature_positions(creature_positions_),
	  creature_fitnesses(std::vector<double>(creature_positions.size(), std::numeric_limits<double>::infinity())),
	  lower_bound(lower_bound_),
	  upper_bound(upper_bound_),
	  survival_rate(survival_rate_),
	  mutation_rate(mutation_rate_),
	  best_creature_index(0),
	  func(func_) {
	assert(world_rank >= 0);
	assert(creature_positions_.size() > 0);
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
	std::cout << "[" << world_rank << "] TODO: implement sorting" << std::endl;
}

void DistributedGeneticAlgorithm::applyCrossover(const size_t) {
	std::cout << "[" << world_rank << "] TODO: implement crossover" << std::endl;
}

void DistributedGeneticAlgorithm::applyMutation(const size_t) {
	std::cout << "[" << world_rank << "] TODO: implement mutation" << std::endl;
}
