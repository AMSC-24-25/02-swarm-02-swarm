#include "FireflyAlgorithm.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip>

// Constructor
FireflyAlgorithm::FireflyAlgorithm(int nFireflies, int dim, double a, double b, double g, double lower, double upper, size_t s)
	: numFireflies(nFireflies), dimensions(dim), alpha(a), beta(b), gamma(g),lower_bound(lower), upper_bound(upper), seed(s) {
	initializeFireflies();
}


// Set the objective function
void FireflyAlgorithm::setObjectiveFunction(std::function<double(const std::vector<double>&)> func) {
    objectiveFunction = func;
}

// Initialize fireflies
void FireflyAlgorithm::initializeFireflies() {
    fireflies.clear();
    fireflies.resize(numFireflies, Firefly(dimensions));  // Initial fill

	std::uniform_real_distribution<double> dist(lower_bound, upper_bound);

#pragma omp parallel for
    for (int i = 0; i < numFireflies; ++i) {
    	std::mt19937 thread_rng(seed + i); // SEED UNICO PER OGNI FIRELY
    	std::vector<double> pos(dimensions);
    	for (int d = 0; d < dimensions; ++d) {
    		pos[d] = dist(thread_rng);
    	}
        fireflies[i] = Firefly(dimensions);
    	fireflies[i].setPosition(pos);
        fireflies[i].setBrightness(std::numeric_limits<double>::max());
    }
}

// Compute Euclidean distance between two fireflies
double FireflyAlgorithm::euclideanDistance(const Firefly& a, const Firefly& b) {
    std::vector<double> posA = a.getPosition();
    std::vector<double> posB = b.getPosition();

    double sum = 0.0;
    for (size_t i = 0; i < posA.size(); ++i) {
        sum += (posA[i] - posB[i]) * (posA[i] - posB[i]);
    }
    return std::sqrt(sum);
}

// Update fireflies based on attractiveness and random perturbation
void FireflyAlgorithm::updateFireflies() {
    std::random_device rd;
    std::mt19937 gen(rd()); //per riproducibilita'con solo seed
    std::uniform_real_distribution<> randomNoise(-1.0, 1.0);

    for (int i = 0; i < numFireflies; ++i) {
        for (int j = 0; j < numFireflies; ++j) {
            if (fireflies[j].getBrightness() < fireflies[i].getBrightness()) {
                double r = euclideanDistance(fireflies[i], fireflies[j]);
                double beta_effective = beta * exp(-gamma * r * r);

                std::vector<double> newPosition = fireflies[i].getPosition();
                std::vector<double> targetPosition = fireflies[j].getPosition();

                for (int d = 0; d < dimensions; ++d) {
                    newPosition[d] += beta_effective * (targetPosition[d] - newPosition[d])
                                    + alpha * randomNoise(gen);
                }
            	// CLIPPING
            	for (int d = 0; d < dimensions; ++d) {
            		newPosition[d] = std::max(lower_bound, std::min(upper_bound, newPosition[d]));
            	}

                fireflies[i].setPosition(newPosition);
            }
        }
    }
}

// Main optimization loop
std::vector<double> FireflyAlgorithm::optimize(int maxIterations) {
	auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < maxIterations; ++iter) {

#pragma omp parallel for
        for (Firefly& firefly : fireflies) {
            firefly.setBrightness(objectiveFunction(firefly.getPosition()));
        }

        updateFireflies();


    	//Trova la firefly migliore dell'iterazione
    	const Firefly* best = &fireflies[0];
    	for (const auto& f : fireflies) {
    		if (f.getBrightness() < best->getBrightness()) {
    			best = &f;
    		}
    	}

    	//Output formattato
    	std::cout << "Iteration n. " << iter + 1 << " / " << maxIterations << "\n";
    	std::cout << "  Current minimum:\n";
    	std::cout << "  f(";
    	const std::vector<double>& pos = best->getPosition();
    	for (size_t i = 0; i < pos.size(); ++i) {
    		std::cout << std::scientific << std::setprecision(6) << pos[i];
    		if (i < pos.size() - 1) std::cout << ", ";
    	}
    	std::cout << ") = " << std::scientific << std::setprecision(6)
				  << best->getBrightness() << "\n";


    }

    // Return the best firefly found
    int bestIndex = 0;
    for (int i = 1; i < numFireflies; ++i) {
        if (fireflies[i].getBrightness() < fireflies[bestIndex].getBrightness()) {
            bestIndex = i;
        }
    }
	auto end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(end - start).count();
	std::cout << "Total execution time: " << std::fixed << std::setprecision(6)
			  << elapsed << " seconds" << std::endl;
	std::cout << std::endl;
    return fireflies[bestIndex].getPosition();
}
