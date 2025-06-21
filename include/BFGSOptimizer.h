#ifdef USE_EIGEN
#ifndef BFGS_OPTIMIZER_H
#define BFGS_OPTIMIZER_H

#include <Eigen/Dense>
#include <functional>
#include <vector>

class BFGSOptimizer {
public:
    BFGSOptimizer(int dimensions);
    void setObjective(std::function<double(const std::vector<double>&)> f);
    std::vector<double> optimize(const std::vector<double>& x0, int maxIterations = 100);

private:
    std::function<double(const std::vector<double>&)> objective;
    int dim;

    Eigen::VectorXd numericalGradient(const Eigen::VectorXd& x);
    double evaluate(const Eigen::VectorXd& x);
};

#endif // BFGS_OPTIMIZER_H
#endif // USE_EIGEN
