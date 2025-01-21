#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>

#include "Particle.hpp"

Particle::Particle(const int dimensions_, const double lower_bound, const double upper_bound, const size_t seed) {
	assert(dimensions_ > 0);

	dimensions = dimensions_;
	bestFitness = std::numeric_limits<double>::infinity();
	std::vector<double> position_(dimensions);
	std::vector<double> velocity_(dimensions);
	std::vector<double> bestLocalPosition_(dimensions);

	std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> position_dist{lower_bound, upper_bound};
	std::uniform_real_distribution<double> velocity_dist{-(upper_bound - lower_bound) * 0.1,
														 (upper_bound - lower_bound) * 0.1};

	for (int i = 0; i < dimensions; i++) {
		position_[i] = position_dist(rnd);
		velocity_[i] = velocity_dist(rnd);
		bestLocalPosition_[i] = position_[i];
	}

	position = position_;
	velocity = velocity_;
	bestLocalPosition = bestLocalPosition_;
}

void Particle::update(const ObjectiveFunction& func, const std::vector<double>& globalBestPosition,
					  const double lower_bound, const double upper_bound, const double c1, const double c2,
					  const double w, const size_t seed) {
	assert(w > 0.0);
	assert(w <= 1.0);

	std::mt19937 rnd{seed};

	std::uniform_real_distribution<double> r{0.0, 1.0};

	for (int i = 0; i < dimensions; i++) {
		const double r1 = r(rnd);
		const double r2 = r(rnd);

		velocity[i] = velocity[i] * w + c1 * r1 * (bestLocalPosition[i] - position[i]) +
					  c2 * r2 * (globalBestPosition[i] - position[i]);

		// solid "sticky" bounds, when a particle reaches a boundary, it sticks
		// to it
		position[i] = std::clamp(position[i] + velocity[i], lower_bound, upper_bound);
	}

	const double newVal = func(position);
	if (newVal < bestFitness) {
		bestFitness = newVal;
		for (int i = 0; i < dimensions; i++) {
			bestLocalPosition[i] = position[i];
		}
	}
}
