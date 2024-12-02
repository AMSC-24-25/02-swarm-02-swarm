#include "Swarm.hpp"


    Swarm::Swarm(vector<Particle>& particles_, vector<double> lower_, vector<double> upper_) {
      particles = particles_;
      minimus = particles[0].bestFitness;
      bestGlobalPosition = particles[0].bestLocalPosition;
      lower = lower_;
      upper = upper_;
    }

    double Swarm::findbestFitness() {
        for (int i = 0; i < particles.size(); ++i) {


            if(particles[i].bestFitness < minimus) {
              minimus = particles[i].bestFitness;
              bestGlobalPosition= particles[i].position;
            }
        }
        return minimus;
    }

    void Swarm::updateParticles() {

      for (int i = 0; i < particles.size(); ++i) {
     
        particles[i].update(bestGlobalPosition);
      }
    }
