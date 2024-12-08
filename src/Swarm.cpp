#include <vector>
#include <omp.h>

#include "Swarm.hpp"

Swarm::Swarm(std::vector<Particle>& particles_, std::vector<double> lower_, std::vector<double> upper_ , double c1_, double c2_, double w_) { 
	particles = particles_;
	minimum = particles[0].bestFitness;
	bestGlobalPosition = particles[0].bestLocalPosition;
	lower = lower_;
	upper = upper_;
	c1 = c1_;
	c2 = c2_;
	w = w_;
}

double Swarm::findbestFitness() {
	double local_minimum = minimum;
	std::vector<double> local_bestPosition = bestGlobalPosition;

#pragma omp parallel
	{	//for avoid race-condition
		double thread_minimum = local_minimum;
		std::vector<double> thread_bestPosition = local_bestPosition;

		//no need to wait all threads
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
#pragma omp parallel for
	for (size_t i = 0; i < particles.size(); ++i) {
		particles[i].update(bestGlobalPosition, lower, upper, c1, c2, w);
	}
}

void Swarm::updateInertia(int max_iterations, double w_min, double w_max){
	w = w - ((w_max - w_min)/ max_iterations );
}
