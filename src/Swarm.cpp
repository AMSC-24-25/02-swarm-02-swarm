#include <vector>
#include <cassert>
#include <omp.h>

#include "Swarm.hpp"

Swarm::Swarm(const std::vector<Particle>& particles_, const double lower_bound_, const double upper_bound_,
			 const double c1_, const double c2_, const double w_, const size_t seed_, ObjectiveFunction& func_,
			 const size_t n_threads_)
	: particles(particles_),
	  bestGlobalPosition(particles[0].bestLocalPosition),
	  minimum(particles[0].bestFitness),
	  lower_bound(lower_bound_),
	  upper_bound(upper_bound_),
	  c1(c1_),
	  c2(c2_),
	  w(w_),
	  seed(seed_),
	  n_threads(n_threads_),
	  func(func_) {}

void Swarm::findBestFitness() {
	double local_minimum = minimum;
	std::vector<double> local_bestPosition = bestGlobalPosition;

#pragma omp parallel num_threads(n_threads) default(none) shared(local_minimum, local_bestPosition)
	{  // for avoid race-condition
		double thread_minimum = local_minimum;
		std::vector<double> thread_bestPosition = local_bestPosition;

		// no need to wait all threads
#pragma omp for nowait
		for (size_t i = 0; i < particles.size(); i++) {
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
	std::copy(local_bestPosition.begin(), local_bestPosition.end(), bestGlobalPosition.begin());
	// bestGlobalPosition = local_bestPosition;
}

void Swarm::updateParticles() {
#pragma omp parallel for schedule(static) num_threads(n_threads)
	for (size_t i = 0; i < particles.size(); i++) {
		// Each particle receives a unique seed (different from the global one) so that each has a different sequence of
		// random numbers
		particles[i].update(func, bestGlobalPosition, lower_bound, upper_bound, c1, c2, w, seed + i);
	}

	// Update seed
	seed += particles.size();
}

void Swarm::updateInertia(const size_t max_iterations, const double w_min, const double w_max) {
	assert(max_iterations > 0);
	assert(w_min > 0.0);
	assert(w_min < w_max);
	assert(w_max <= 1.0);
	w = w - ((w_max - w_min) / static_cast<double>(max_iterations));
}
