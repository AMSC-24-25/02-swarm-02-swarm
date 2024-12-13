#include <vector>
#include <cassert>
#include <omp.h>

#include "Swarm.hpp"

Swarm::Swarm(const std::vector<Particle>& particles_, const std::vector<double>& lower_,
			 const std::vector<double>& upper_, const double c1_, const double c2_, const double w_, const size_t seed_,
			 ObjectiveFunction& func_)
	: particles(particles_),
	  bestGlobalPosition(particles[0].bestLocalPosition),
	  lower(lower_),
	  upper(upper_),
	  minimum(particles[0].bestFitness),
	  c1(c1_),
	  c2(c2_),
	  w(w_),
	  seed(seed_),
	  func(func_) {}

double Swarm::findbestFitness() {
	double local_minimum = minimum;
	std::vector<double> local_bestPosition = bestGlobalPosition;

#pragma omp parallel
	{  // for avoid race-condition
		double thread_minimum = local_minimum;
		std::vector<double> thread_bestPosition = local_bestPosition;

		// no need to wait all threads
#pragma omp for nowait
		for (size_t i = 0; i < particles.size(); ++i) {
			if (particles[i].bestFitness < thread_minimum) {
				thread_minimum = particles[i].bestFitness;
				thread_bestPosition = particles[i].position;
			}
		}

#pragma omp critical
		{
			if (thread_minimum < local_minimum) {
				local_minimum = thread_minimum;
				local_bestPosition = thread_bestPosition;
			}
		}
	}

	minimum = local_minimum;
	bestGlobalPosition = local_bestPosition;
	return minimum;
}

void Swarm::updateParticles() {
#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < particles.size(); ++i) {
		// Each particle receives a unique seed (different from the global one) so that each has a different sequence of
		// random numbers
		particles[i].update(func, bestGlobalPosition, lower, upper, c1, c2, w, seed + i + 1);
	}
}

void Swarm::updateInertia(const int max_iterations, const double w_min, const double w_max) {
	assert(max_iterations > 0);
	assert(w_min > 0.0);
	assert(w_min < w_max);
	assert(w_max <= 1.0);
	w = w - ((w_max - w_min) / max_iterations);
}
