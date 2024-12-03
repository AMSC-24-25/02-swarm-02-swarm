#include <vector>
#include <iostream>
#include <iomanip>

#include "Swarm.hpp"

int main() {
	const int dimensions = 2;
	const int num_particles = 100;
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions);
	std::vector<double> upperBound(dimensions);

	lowerBound[0] = -100.0;
	lowerBound[1] = -100.0;
	upperBound[0] = 100.0;
	upperBound[1] = 100.0;

	for (int i = 0; i < num_particles; ++i) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound));
	}

	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound);

	const int max_iterations = 100;

	for (int i = 0; i < max_iterations; ++i) {
		swarm.updateParticles();
		swarm.findbestFitness();
		std::cout << "Iteration n. " << i << " / " << max_iterations
				  << std::endl;
		std::cout << "  Current minimum: " << std::scientific << swarm.minimum
				  << std::endl;
		std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Minimum found:" << std::endl;
	std::cout << "  f(";
	for (int i = 0; i < dimensions; i++) {
		std::cout << std::scientific << swarm.bestGlobalPosition[i];
		if (i < dimensions - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << std::scientific << swarm.minimum << std::endl;
	std::cout << std::endl;

	return 0;
}
