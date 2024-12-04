#include <vector>
#include <limits>
#include <random>
#include <algorithm>

#include "Particle.hpp"

Particle::Particle(int dimensions_, std::vector<double>& lower, std::vector<double>& upper) {
	dimensions = dimensions_;
	bestFitness = std::numeric_limits<double>::infinity();
	std::vector<double> position_(dimensions);
	std::vector<double> velocity_(dimensions);
	std::vector<double> bestLocalPosition_(dimensions);

	std::random_device dev;
	std::mt19937 rnd{dev()};

	for (int i = 0; i < dimensions; ++i) {
		std::uniform_real_distribution<double> position_dist{lower[i], upper[i]};
		position_[i] = position_dist(rnd);
		std::uniform_real_distribution<double> velocity_dist{-(upper[i] - lower[i]) * 0.1, (upper[i] - lower[i]) * 0.1};
		velocity_[i] = velocity_dist(rnd);
		bestLocalPosition_[i] = position_[i];
	}

	position = position_;
	velocity = velocity_;
	bestLocalPosition = bestLocalPosition_;
}

void Particle::update(std::vector<double>& globalBestPosition, const std::vector<double>& lower,
					  const std::vector<double>& upper, double c1, double c2, double w) {

	double r1;
	double r2;
	std::random_device dev;
	std::mt19937 rnd{dev()};

	std::uniform_real_distribution<double> r{0, 1};

	for (int i = 0; i < dimensions; ++i) {

		
		r1 = r(rnd);
		r2= r(rnd);

		velocity[i] = velocity[i] * w + c1 * r1 * (bestLocalPosition[i] - position[i]) +
					  c2 * r2 * (globalBestPosition[i] - position[i]);

		// solid "sticky" bounds, when a particle reaches a boundary, it sticks
		// to it
		position[i] = std::clamp(position[i] + velocity[i], lower[i], upper[i]);
	}
	const double newVal = ObjectiveFunction().getValueFunction(position);
	if (newVal < bestFitness) {
		bestFitness = newVal;

		for (int i = 0; i < dimensions; ++i) {
			bestLocalPosition[i] = position[i];
		}
	}
}
