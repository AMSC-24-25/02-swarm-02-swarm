#include <vector>
#include <limits>

#include "Particle.hpp"

Particle::Particle(int dimensions_, std::vector<double>& lower,
				   std::vector<double>& upper) {
	dimensions = dimensions_;
	bestFitness = std::numeric_limits<double>::infinity();
	std::vector<double> position_(dimensions);
	std::vector<double> velocity_(dimensions);
	std::vector<double> bestLocalPosition_(dimensions);

	for (int i = 0; i < dimensions; ++i) {
		position_[i] = (lower[i] + static_cast<double>(rand()) / RAND_MAX *
									   (upper[i] - lower[i]));
		velocity_[i] = ((static_cast<double>(rand()) / RAND_MAX * 2 - 1) *
						(upper[i] - lower[i]) * 0.1);
		bestLocalPosition_[i] = (position_[i]);
	}
	position = position_;
	velocity = velocity_;
	bestLocalPosition = bestLocalPosition_;
}

void Particle::update(std::vector<double>& globalBestPosition) {
	for (int i = 0; i < dimensions; ++i) {
		velocity[i] = velocity[i] * (0.4) +
					  (0.6) * (bestLocalPosition[i] - position[i]) +
					  (0.5) * (globalBestPosition[i] - position[i]);
		// sistemare quando esce da lower/upper bound
		position[i] = position[i] + velocity[i];
	}
	double newValueOfFunction = ObjectiveFunction().getValueFunction(position);
	if (newValueOfFunction < bestFitness) {
		bestFitness = newValueOfFunction;

		for (int i = 0; i < dimensions; ++i) {
			bestLocalPosition[i] = position[i];
		}
	}
}
