//
// Created by javed-abdullah on 12/2/24.
//

#include "Swarm.hpp"
#define DIM 2

int main() {

    std::vector<Particle>  swarmParticles;
    std::vector<double> lowerBound;
    std::vector<double> upperBound;

    lowerBound.push_back(-100);
    lowerBound.push_back(-100);
    upperBound.push_back(100);
    upperBound.push_back(100);

    std::cout<<lowerBound.size()<<std::endl;
    std::cout<<upperBound.size()<<std::endl;

    for (int i = 0; i < 100; ++i) {
        swarmParticles.push_back({DIM,lowerBound,upperBound});
    }
    Swarm swarm = Swarm(swarmParticles);

    for (int i = 0; i < 100; ++i) {
      swarm.updateParticles();
       swarm.findbestFitness();
       cout << swarm.minimus << endl;
    }




  return 0;
}