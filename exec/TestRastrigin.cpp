#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "Swarm.hpp"
#include "Rastrigin.hpp"

double absolute_error(const double expected, const double actual) { 
    return std::abs(expected - actual); 
}

double calculate_max_distance(const std::vector<double>& lowerBound, const std::vector<double>& upperBound) {
    // Calcola la distanza massima possibile dal minimo globale
    double max_distance = 0.0;
    for (size_t i = 0; i < lowerBound.size(); ++i) {
        max_distance += std::pow(std::max(std::abs(lowerBound[i]), std::abs(upperBound[i])), 2);
    }
    return std::sqrt(max_distance);
}

int main() {
    const int dimensions = 6;
    const int num_particles = 10000;
    const int max_iterations = 100;
    double distance_from_globalminimum = 0.0;

    std::vector<Particle> swarmParticles;
    std::vector<double> lowerBound(dimensions, -5.12);
    std::vector<double> upperBound(dimensions, 5.12);

    // Calcola la distanza massima possibile dal minimo globale
    const double max_distance = calculate_max_distance(lowerBound, upperBound);

    // Inizializza il primo sciame di particelle
    for (int i = 0; i < num_particles; ++i) {
        swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, 42));
    }

    const double w_max = 0.9;
    const double w_min = 0.4;
    const double w = w_max;

    Rastrigin r;
    Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, 2.0, 2.0, w, 42, r, 1);

    const auto beginning = omp_get_wtime();

    // Esegui l'ottimizzazione per il primo tentativo
    for (int i = 0; i < max_iterations; ++i) {
        swarm.updateInertia(max_iterations, w_min, w_max);
        swarm.updateParticles();
        swarm.findBestFitness();
    }
    const auto end = omp_get_wtime();
    std::cout << "Sequential attempt: " << std::endl;
    std::cout << "Residual: " << swarm.minimum - 0.0 << std::endl;
    std::cout << "Total execution time: " << (end - beginning) << " seconds" << std::endl;

    // Calcola la distanza dal minimo globale
    for (int i = 0; i < dimensions; i++) {
        distance_from_globalminimum += (swarm.bestGlobalPosition[i]) * (swarm.bestGlobalPosition[i]);
    }
    double distance = sqrt(distance_from_globalminimum);

    // Calcola la distanza normalizzata
    double normalized_distance = distance / max_distance;
    std::cout << "Distance from the global minimum: " << distance << std::endl;
    std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;

    std::cout << "******************************************************" << std::endl;


    std::vector<Particle> swarmParticles2;

    for (int i = 0; i < num_particles; ++i) {
        swarmParticles2.push_back(Particle(dimensions, lowerBound, upperBound, 42));
    }

    Swarm swarm2 = Swarm(swarmParticles2, lowerBound, upperBound, 2.0, 2.0, w, 42, r, 100);

    const auto beginning2 = omp_get_wtime();
    for (int i = 0; i < max_iterations; ++i) {
        swarm2.updateInertia(max_iterations, w_min, w_max);
        swarm2.updateParticles();
        swarm2.findBestFitness();
    }
    const auto end2 = omp_get_wtime();

    std::cout << "Parallel attempt: " << std::endl;
    std::cout << "Residual: " << swarm2.minimum - 0.0 << std::endl;
    std::cout << "Total execution time: " << (end2 - beginning2) << " seconds" << std::endl;

    distance_from_globalminimum = 0.0;
    for (int i = 0; i < dimensions; i++) {
        distance_from_globalminimum += (swarm2.bestGlobalPosition[i]) * (swarm2.bestGlobalPosition[i]);
    }
    distance = sqrt(distance_from_globalminimum);

    // Calcola la distanza normalizzata
    normalized_distance = distance / max_distance;
    std::cout << "Distance from the global minimum: " << distance << std::endl;
    std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;

    std::cout << "*****************************************************" << std::endl;

    distance_from_globalminimum = 0.0;

    std::vector<Particle> swarmParticles3;

    for (int i = 0; i < num_particles; ++i) {
        swarmParticles3.push_back(Particle(dimensions, lowerBound, upperBound, 42));
    }

    const double w_max2 = 0.8;
    const double w_min2 = 0.4;
    const double w2 = w_max2;

    Swarm swarm3 = Swarm(swarmParticles3, lowerBound, upperBound, 2, 2, w2, 42, r, 100);

    const auto beginning3 = omp_get_wtime();

    for (int i = 0; i < max_iterations; ++i) {
        swarm3.updateInertia(max_iterations, w_min2, w_max2);
        swarm3.updateParticles();
        swarm3.findBestFitness();
    }
    const auto end3 = omp_get_wtime();
    std::cout << "Best attempt (parallel): " << std::endl;
    std::cout << "Residual: " << swarm3.minimum - 0.0 << std::endl;
    std::cout << "Total execution time: " << (end3 - beginning3) << " seconds" << std::endl;

    for (int i = 0; i < dimensions; i++) {
        distance_from_globalminimum += (swarm3.bestGlobalPosition[i]) * (swarm3.bestGlobalPosition[i]);
    }
    distance = sqrt(distance_from_globalminimum);

    // Calcola la distanza normalizzata
    normalized_distance = distance / max_distance;
    std::cout << "Distance from the global minimum: " << distance << std::endl;
    std::cout << "Normalized distance (percentage): " << normalized_distance * 100 << "%" << std::endl;

    std::cout << "In the best attempt we used: w_max = " << w_max2 << ", w_min = " << w_min2 << ", c1 = " << 2 << ", c2 = " << 2 << std::endl;

    return 0;
}
