//
// Created by javed-abdullah on 12/2/24.
//
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>


using namespace std;


// Funzione obiettivo (esempio: funzione Sphere) (sommatoria delle x^2)
double objectiveFunction(const vector<double>& position) {
    double sum = 0.0;
    for (double x : position) {
        sum += x * x;
    }
    return sum;
}


// Classe per rappresentare una particella
class Particle {
public:
    vector<double> position;    // Posizione attuale
    vector<double> velocity;    // Velocità attuale
    vector<double> bestLocalPosition; // Miglior posizione personale
    double bestFitness;         // Miglior fitness personale (f(bestLocalpostion)
    int dimensions, lower, upper;

    Particle(int dimensions, double lowerBound, double upperBound) {
      this->dimensions = dimensions;
      this->lower= lowerBound;
       this->upper = upperBound;
        position.resize(dimensions);
        velocity.resize(dimensions);
       bestLocalPosition.resize(dimensions);
        bestFitness = numeric_limits<double>::infinity();

        for (int i = 0; i < dimensions; ++i) {
            position[i] = lowerBound + static_cast<double>(rand()) / RAND_MAX * (upperBound - lowerBound);

            velocity[i] = (static_cast<double>(rand()) / RAND_MAX * 2 - 1) * (upperBound - lowerBound) * 0.1;
            bestLocalPosition[i] = position[i];
        }
    }
    void update(vector<double>& globalBestPosition) {
        for (int i = 0; i < dimensions; ++i) {
            velocity[i] =  velocity[i]*(0.4) + (0.6)*(bestLocalPosition[i]- position[i]) + (0.5) * (globalBestPosition[i] - position[i]);
            //sistemrae quando esce da lower/upper bound
            position[i] = position[i] + velocity[i];

        }
        double newValuoOfFunction = objectiveFunction(position);
        if(newValuoOfFunction < bestFitness) {
          bestFitness = newValuoOfFunction;

          for(int i = 0; i < dimensions; ++i) {
             bestLocalPosition[i] = position[i];
          }
        }

    }


};

class Swarm{

  public:
    vector<Particle> particles;
    vector<double> bestGlobalPosition;
    double minimus;
    Swarm(vector<Particle>& particles) {
       this-> particles = particles;
      this->  minimus = particles[0].bestFitness;
      bestGlobalPosition = particles[0].bestLocalPosition;
    }

    double findbestFitness() {
        for (int i = 0; i < particles.size(); ++i) {


            if(particles[i].bestFitness < minimus) {
              minimus = particles[i].bestFitness;
              bestGlobalPosition= particles[i].position;
            }
        }
        return minimus;
    }

    void updateParticles() {

      for (int i = 0; i < particles.size(); ++i) {
     
        particles[i].update(bestGlobalPosition);
      }
    }
};


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