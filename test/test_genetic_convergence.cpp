#include <vector>
#include <random>

#include <gtest/gtest.h>

#include "GeneticAlgorithm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) {
	assert(std::isfinite(expected));
	assert(std::isfinite(actual));
	return std::abs(expected - actual);
}

TEST(GeneticConvergence, Sphere) {
	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	std::vector<Creature> creatures;

	std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
	for (size_t i{0}; i < num_creatures; i++) {
		std::vector<double> tmp(dimensions);
		std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
		creatures.push_back(Creature(tmp));
	}

	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	Sphere s;

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, s, 1);

	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();
	}

	EXPECT_LE(absolute_error(0.0, ga.bestCreature.fitness), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), ga.bestCreature.position.at(i)), 1e-3);
	}
}

TEST(GeneticConvergence, EuclideanDistance) {
	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	std::vector<Creature> creatures;

	std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
	for (size_t i{0}; i < num_creatures; i++) {
		std::vector<double> tmp(dimensions);
		std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
		creatures.push_back(Creature(tmp));
	}

	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	EuclideanDistance ed;

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, ed, 1);

	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();
	}

	EXPECT_LE(absolute_error(0.0, ga.bestCreature.fitness), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), ga.bestCreature.position.at(i)), 1e-3);
	}
}

TEST(GeneticConvergence, Rosenbrock) {
	const size_t dimensions = 2;
	const size_t num_creatures = 1'000;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	std::vector<Creature> creatures;

	std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
	for (size_t i{0}; i < num_creatures; i++) {
		std::vector<double> tmp(dimensions);
		std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
		creatures.push_back(Creature(tmp));
	}

	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	Rosenbrock r;

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, r, 1);

	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();
	}

	EXPECT_LE(absolute_error(0.0, ga.bestCreature.fitness), 1e-3);

	std::vector<double> expected_minimum(dimensions);
	for (size_t i = 0; i < dimensions; i++) {
		expected_minimum[i] = std::pow(r.a, static_cast<double>(i + 1));
	}
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), ga.bestCreature.position.at(i)), 0.1);
	}
}

TEST(GeneticConvergence, Rastrigin) {
	const size_t dimensions = 2;
	const size_t num_creatures = 100;
	const size_t max_iterations = 1'000;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	std::vector<Creature> creatures;

	std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
	for (size_t i{0}; i < num_creatures; i++) {
		std::vector<double> tmp(dimensions);
		std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
		creatures.push_back(Creature(tmp));
	}

	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	Rastrigin r;

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, r, 1);

	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();
	}

	EXPECT_LE(absolute_error(0.0, ga.bestCreature.fitness), 1e-3);

	std::vector<double> expected_minimum(dimensions, 0.0);
	for (size_t i = 0; i < dimensions; i++) {
		EXPECT_LE(absolute_error(expected_minimum.at(i), ga.bestCreature.position.at(i)), 1e-3);
	}
}
