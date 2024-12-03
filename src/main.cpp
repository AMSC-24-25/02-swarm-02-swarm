#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "Swarm.hpp"

void print_minimum(const Swarm& swarm, const size_t dimensions) {
	std::cout << "  f(";
	for (size_t i = 0; i < dimensions; i++) {
		std::cout << std::scientific << swarm.bestGlobalPosition[i];
		if (i < dimensions - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << std::scientific << swarm.minimum << std::endl;
}

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
	const auto beginning = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < max_iterations; ++i) {
		const auto start_iteration = std::chrono::high_resolution_clock::now();
		swarm.updateParticles();
		swarm.findbestFitness();
		const auto end_iteration = std::chrono::high_resolution_clock::now();
		std::cout << "Iteration n. " << i << " / " << max_iterations << std::endl;
		std::cout << "  Current minimum: " << std::endl;
		print_minimum(swarm, dimensions);
		std::cout
			<< "  Execution time: " << std::fixed << std::setprecision(6)
			<< (static_cast<double>(
					std::chrono::duration_cast<std::chrono::nanoseconds>(end_iteration - start_iteration).count()) *
				1e-9)
			<< " seconds" << std::endl;
		std::cout << std::endl;
	}

	const auto end = std::chrono::high_resolution_clock::now();

	std::cout << std::endl;
	std::cout << "Minimum found:" << std::endl;
	print_minimum(swarm, dimensions);
	std::cout << "  Total execution time: " << std::fixed << std::setprecision(6)
			  << (static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - beginning).count()) *
				  1e-9)
			  << " seconds" << std::endl;
	std::cout << std::endl;

	return 0;
}
