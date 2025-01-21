#include <vector>
#include <cassert>
#include <cmath>

#include "GeneticAlgorithm.hpp"

GeneticAlgorithm::GeneticAlgorithm(const std::vector<Creature>& creatures_, const double lower_bound_,
								   const double upper_bound_, const double mutation_rate_, const double crossover_rate_,
								   const double survival_rate_)
	: creatures(creatures_),
	  lower_bound(lower_bound_),
	  upper_bound(upper_bound_),
	  survival_rate(survival_rate_),
	  mutation_rate(mutation_rate_),
	  crossover_rate(crossover_rate_),
	  bestCreature(creatures[0]),
	  bestFitness(std::numeric_limits<double>::infinity()) {
	assert(creatures.size() > 0);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);
	assert(survival_rate > 0.0 && survival_rate < 1.0);
	assert(mutation_rate >= 0.0 && mutation_rate < 1.0);
	assert(crossover_rate >= 0.0 && crossover_rate < 1.0);
}
