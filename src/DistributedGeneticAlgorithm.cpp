#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>

#include "DistributedGeneticAlgorithm.hpp"

DistributedGeneticAlgorithm::DistributedGeneticAlgorithm(const std::vector<Creature>& creatures_,
														 const double lower_bound_, const double upper_bound_,
														 const double survival_rate_, const double mutation_rate_,
														 ObjectiveFunction& func_)
	: creatures(creatures_),
	  lower_bound(lower_bound_),
	  upper_bound(upper_bound_),
	  survival_rate(survival_rate_),
	  mutation_rate(mutation_rate_),
	  bestCreature(creatures[0]),
	  func(func_) {
	assert(creatures.size() > 0);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);
	assert(survival_rate > 0.0 && survival_rate < 1.0);
	assert(mutation_rate >= 0.0 && mutation_rate < 1.0);
}

void DistributedGeneticAlgorithm::evaluateCreatures() {
	std::cout << "TODO: implement evaluation" << std::endl;
}

void DistributedGeneticAlgorithm::sortCreatures() {
	std::cout << "TODO: implement sorting" << std::endl;
}

void DistributedGeneticAlgorithm::applyCrossover(const size_t) {
	std::cout << "TODO: implement crossover" << std::endl;
}

void DistributedGeneticAlgorithm::applyMutation(const size_t) {
	std::cout << "TODO: implement mutation" << std::endl;
}
