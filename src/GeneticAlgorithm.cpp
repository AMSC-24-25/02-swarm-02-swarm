#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <limits>
#include <utility>
#include <numeric>

#include "GeneticAlgorithm.hpp"

GeneticAlgorithm::GeneticAlgorithm(const std::vector<Creature>& creatures_, const double lower_bound_,
								   const double upper_bound_, const double survival_rate_, const double mutation_rate_,
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

void GeneticAlgorithm::evaluateCreatures() {
	for (size_t i{0}; i < creatures.size(); i++) {
		// Avoid evaluating the same creatures multiple times
		if (std::isfinite(creatures.at(i).fitness)) {
			continue;
		}

		creatures.at(i).fitness = func(creatures.at(i).position);
	}
}

void GeneticAlgorithm::sortCreatures() {
	std::sort(creatures.begin(), creatures.end(),
			  [](const Creature& a, const Creature& b) { return a.fitness > b.fitness; });

	bestCreature = creatures.at(0);
}

void GeneticAlgorithm::applyCrossover() {
	const size_t n_creatures = creatures.size();
	const size_t survived = creatures.size() * survival_rate;

	// The population which did not survive is not actually removed, it gets overwritten during crossover

	// Fill the rest of the population with new offspring
	std::mt19937 rnd{42};
	std::uniform_int_distribution<size_t> index_dist{0, survived - 1};
	std::bernoulli_distribution bool_dist{0.5};
	for (size_t i{survived}; i < n_creatures; i++) {
		// Choosing father randomly
		const size_t father_index = index_dist(rnd);
		assert(father_index < survived);
		const Creature& father = creatures.at(father_index);

		// Choose mother randomly
		size_t mother_index;
		do {
			mother_index = index_dist(rnd);
		} while (father_index == mother_index);
		assert(mother_index < survived);
		assert(father_index != mother_index);
		const Creature& mother = creatures.at(mother_index);

		for (size_t j{0}; j < creatures.at(i).position.size(); j++) {
			creatures.at(i).position.at(j) = bool_dist(rnd) ? father.position.at(j) : mother.position.at(j);
		}
		creatures.at(i).fitness = std::numeric_limits<double>::infinity();
	}

	// Make sure that no creature is added nor removed from the vector during crossover
	assert(creatures.size() == n_creatures);
}

void GeneticAlgorithm::applyMutation() {
	std::mt19937 rnd{42};
	std::bernoulli_distribution bool_dist{mutation_rate};
	std::uniform_real_distribution<double> position_dist{lower_bound, upper_bound};

	const size_t dimensions = creatures.at(0).position.size();
	std::uniform_int_distribution<size_t> index_dist{0, dimensions - 1};

	for (Creature& c : creatures) {
		if (!bool_dist(rnd)) {
			continue;
		}

		c.position[index_dist(rnd)] = position_dist(rnd);
	}
}
