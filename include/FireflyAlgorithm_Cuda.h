#ifndef FIREFLY_ALGORITHM_CUDA_H
#define FIREFLY_ALGORITHM_CUDA_H

#include "FireflyAlgorithm.h" // NON solo Firefly, ereditiamo da questa
#include <vector>
#include <functional>

// 🔹 Wrapper CUDA kernel
extern "C" void launchUpdateFirefliesCUDA(
    double* d_positions, double* d_brightness,
    int numFireflies, int dimensions,
    double alpha, double beta, double gamma,
    unsigned int seed,
    double lower_bound, double upper_bound,
    int threadsPerBlock
);

enum ObjectiveType {
    SPHERE = 0,
    RASTRIGIN = 1,
    ROSENBROCK = 2
};
enum FitnessMode {
    FITNESS_GPU, //  fitness on GPU
    FITNESS_CPU  //  fitness on CPU
};

// 🔹 Ereditiamo da FireflyAlgorithm
class FireflyAlgorithm_Cuda : public FireflyAlgorithm {
public:
    FireflyAlgorithm_Cuda(int numFireflies, int dimensions, double alpha, double beta, double gamma,  double lower_bound = -5, double upper_bound = 5, size_t seed = 42);

    void setObjectiveFunction(ObjectiveType type); // overload
    void setObjectiveFunction(std::function<double(const std::vector<double>&)> func);

    std::vector<double> optimize(int maxIterations) override;

private:
    ObjectiveType objectiveType = SPHERE;
    FitnessMode fitnessMode = FITNESS_GPU;
    std::function<double(const std::vector<double>&)> objectiveFunctionCPU; // function host-side

};

#endif // FIREFLY_ALGORITHM_CUDA_H
