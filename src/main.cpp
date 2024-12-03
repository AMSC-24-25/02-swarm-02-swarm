//
// Created by javed-abdullah on 12/2/24.
//

#include "Swarm.hpp"

int main() {

    std::vector<Particle>  swarmParticles;

    for (int i = 0; i < 100; ++i) {
        swarmParticles.push_back({1,-1000,1000});
    }

    Swarm swarm = Swarm(swarmParticles);

    for (int i = 0; i < 100; ++i) {
      swarm.updateParticles();
       swarm.findbestFitness();
       cout << swarm.minimus << endl;
    }




  return 0;
}