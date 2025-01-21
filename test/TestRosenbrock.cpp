#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "Swarm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"

double absolute_error(const double expected, const double actual) {
	return std::abs(expected - actual);
}

double calculate_max_distance(const int dimensions, const double lower_bound, const double upper_bound) {
	// Calcola la distanza massima possibile dal minimo globale
	double max_distance = 0.0;
	for (int i = 0; i < dimensions; i++) {
		max_distance += std::pow(std::max(std::abs(lower_bound), std::abs(upper_bound)), 2);
	}
	return std::sqrt(max_distance);
}

int main() {
	const int dimensions = 6;
	const int num_particles = 10000;
	const int max_iterations = 100;
	const double lower_bound = -100.0;
	const double upper_bound = 100.0;
	double distance_from_globalminimum = 0.0;
	std::vector<Particle> swarmParticles;

	const double max_distance = calculate_max_distance(dimensions, lower_bound, upper_bound);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, 42));
	}

	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Rosenbrock r;
	Swarm swarm = Swarm(swarmParticles, lower_bound, upper_bound, 2.0, 2.0, w, 42, r, 1);

	const auto beginning = omp_get_wtime();

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
	}
	const auto end = omp_get_wtime();
	std::cout << "Sequential attempt" << std::endl;
	std::cout << "Residual is: " << swarm.minimum - 0.0 << std::endl;
	std::cout << "Total execution time: " << (end - beginning) << " seconds" << std::endl;

	for (int i = 0; i < dimensions; i++) {
		distance_from_globalminimum =
			distance_from_globalminimum + (swarm.bestGlobalPosition[i]) * (swarm.bestGlobalPosition[i]);
	}
	double distance = std::sqrt(distance_from_globalminimum);

	double normalized_distance = distance / max_distance;

	std::cout << "Distance from the global minimum is: " << distance << std::endl;
	std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;
	std::cout << "******************************************************" << std::endl;

	std::vector<Particle> swarmParticles2;

	for (int i = 0; i < num_particles; i++) {
		swarmParticles2.push_back(Particle(dimensions, lower_bound, upper_bound, 42));
	}

	Swarm swarm2 = Swarm(swarmParticles2, lower_bound, upper_bound, 2.0, 2.0, w, 42, r, 100);

	const auto beginning2 = omp_get_wtime();
	for (int i = 0; i < max_iterations; i++) {
		swarm2.updateInertia(max_iterations, w_min, w_max);
		swarm2.updateParticles();
		swarm2.findBestFitness();
	}
	const auto end2 = omp_get_wtime();

	std::cout << "Parallel attempt" << std::endl;
	std::cout << "The residual is: " << swarm2.minimum - 0.0 << std::endl;
	std::cout << "Total execution time: " << (end2 - beginning2) << " seconds" << std::endl;
	distance_from_globalminimum = 0.0;
	for (int i = 0; i < dimensions; i++) {
		distance_from_globalminimum =
			distance_from_globalminimum + (swarm2.bestGlobalPosition[i]) * (swarm2.bestGlobalPosition[i]);
	}
	distance = std::sqrt(distance_from_globalminimum);

	// Calcola la distanza normalizzata
	normalized_distance = distance / max_distance;

	std::cout << "In the parallel attempt the distance from the global minimum is: " << distance << std::endl;
	std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;

	std::cout << "*****************************************************" << std::endl;

	distance_from_globalminimum = 0.0;

	std::vector<Particle> swarmParticles3;

	for (int i = 0; i < num_particles; i++) {
		swarmParticles3.push_back(Particle(dimensions, lower_bound, upper_bound, 42));
	}

	const double w_max2 = 0.9;
	const double w_min2 = 0.3;
	const double w2 = w_max2;

	Swarm swarm3 = Swarm(swarmParticles3, lower_bound, upper_bound, 1, 1, w2, 42, r, 100);

	const auto beginning3 = omp_get_wtime();

	for (int i = 0; i < max_iterations; i++) {
		swarm3.updateInertia(max_iterations, w_min2, w_max2);
		swarm3.updateParticles();
		swarm3.findBestFitness();
	}
	const auto end3 = omp_get_wtime();
	std::cout << "Best attempt(parallel):" << std::endl;
	std::cout << "In the best attempt (parallel)  the residual is: " << swarm3.minimum - 0.0 << std::endl;
	std::cout << "In the best attempt (parallel) the total execution time : " << (end3 - beginning3) << " seconds"
			  << std::endl;

	for (int i = 0; i < dimensions; i++) {
		distance_from_globalminimum =
			distance_from_globalminimum + (swarm3.bestGlobalPosition[i]) * (swarm3.bestGlobalPosition[i]);
	}
	distance = std::sqrt(distance_from_globalminimum);
	normalized_distance = distance / max_distance;
	std::cout << "In the best attempt (parallel) the distance from the global minimum is: " << distance << std::endl;
	std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;
	std::cout << "In the best attempt we used : w_max = " << w_max2 << ", w_min = " << w_min2 << ", c1 = " << 1
			  << ", c2 = " << 1 << std::endl;

	return 0;
}