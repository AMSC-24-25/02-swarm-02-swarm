#ifndef SWARM_HPP
#define SWARM_HPP

#include "Particle.hpp"

class Swarm{

  public:
    vector<Particle> particles;
    vector<double> bestGlobalPosition;
    double minimus;
    Swarm(vector<Particle>& particles);

    double findbestFitness();

    void updateParticles();
};

#endif