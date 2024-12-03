#include <vector>
#include <iostream>

#include "Swarm.hpp"

int main() {
	const int dimension = 2;
	const int num_particles = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimension);
	std::vector<double> upperBound(dimension);

	lowerBound[0] = -100;
	lowerBound[1] = -100;
	upperBound[0] = 100;
	upperBound[1] = 100;

	for (int i = 0; i < num_particles; ++i) {
		swarmParticles.push_back(Particle(dimension, lowerBound, upperBound));
	}

	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound);

	for (int i = 0; i < 100; ++i) {
		swarm.updateParticles();
		swarm.findbestFitness();
		std::cout << swarm.minimus << std::endl;
	}

	return 0;
}
