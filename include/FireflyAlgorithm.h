#ifndef FIREFLY_ALGORITHM_H
#define FIREFLY_ALGORITHM_H

#include "Firefly.h"
#include <vector>
#include <functional>

class FireflyAlgorithm {
public:
    FireflyAlgorithm(int numFireflies, int dimensions, double alpha, double beta, double gamma);

    void setObjectiveFunction(std::function<double(const std::vector<double>&)> func);
   virtual std::vector<double> optimize(int maxIterations);

protected:
    int numFireflies;
    int dimensions;
    double alpha, beta, gamma;

    std::vector<Firefly> fireflies;
    std::function<double(const std::vector<double>&)> objectiveFunction;

    void initializeFireflies();
    void updateFireflies();
    double euclideanDistance(const Firefly& a, const Firefly& b);
};

#endif
