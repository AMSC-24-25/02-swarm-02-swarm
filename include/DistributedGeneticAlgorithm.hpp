#ifndef DISTRIBUTED_GENETIC_ALGORITHM_HPP
#define DISTRIBUTED_GENETIC_ALGORITHM_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class DistributedGeneticAlgorithm {
   public:
	const int world_rank;
	const int world_size;
	const size_t total_creatures;
	std::vector<std::vector<double>> creature_positions;
	std::vector<double> creature_fitnesses;
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	size_t best_creature_index;
	const ObjectiveFunction& func;

	DistributedGeneticAlgorithm(const int world_rank, const int world_size, const size_t total_creatures,
								const std::vector<std::vector<double>>& creature_positions, const double lower_bound,
								const double upper_bound, const double mutation_rate, const double survival_rate,
								ObjectiveFunction& func);

	void evaluateCreatures();

	void sortCreatures();

	void applyCrossover(const size_t seed);

	void applyMutation(const size_t seed);
};

#endif	// DISTRIBUTED_GENETIC_ALGORITHM_HPP
