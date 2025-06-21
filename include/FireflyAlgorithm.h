#ifndef FIREFLY_ALGORITHM_H
#define FIREFLY_ALGORITHM_H

#include "Firefly.h"
#include <vector>
#include <functional>

class FireflyAlgorithm {
public:
    FireflyAlgorithm(int numFireflies, int dimensions, double alpha, double beta, double gamma,  double lower_bound = -5, double upper_bound = 5, size_t seed = 42);

    void setObjectiveFunction(std::function<double(const std::vector<double>&)> func);
   virtual std::vector<double> optimize(int maxIterations);

protected:
    int numFireflies;
    int dimensions;
    double alpha, beta, gamma;
	double lower_bound = -5.0;
	double upper_bound = 5.0;
	size_t seed = 42;

    std::vector<Firefly> fireflies;
    std::function<double(const std::vector<double>&)> objectiveFunction;

    void initializeFireflies();
    void updateFireflies();
    double euclideanDistance(const Firefly& a, const Firefly& b);
};

#endif
