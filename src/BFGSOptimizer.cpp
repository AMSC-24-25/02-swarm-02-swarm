#ifdef USE_EIGEN

#include "BFGSOptimizer.h"
#include <iostream>

BFGSOptimizer::BFGSOptimizer(int dimensions) : dim(dimensions) {}

void BFGSOptimizer::setObjective(std::function<double(const std::vector<double>&)> f) {
    // Move assignment: avoids copying the objective function (e.g. heavy lambda captures)
    // Transfer ownership of the function 'f' to the class member 'objective'
    objective = std::move(f);

}

double BFGSOptimizer::evaluate(const Eigen::VectorXd& x) {
    std::vector<double> x_std(x.data(), x.data() + x.size()); //used for convert the type for objective
    return objective(x_std);
}

Eigen::VectorXd BFGSOptimizer::numericalGradient(const Eigen::VectorXd& x) {
    /*
    * Computes the gradient of the objective function numerically using central differences:
    * df/dx_i ≈ [f(x + εe_i) - f(x - εe_i)] / (2ε)
     */


    const double eps = 1e-6;
    Eigen::VectorXd grad(dim);
    for (int i = 0; i < dim; ++i) {
        Eigen::VectorXd x1 = x, x2 = x;  // Create two perturbed versions of x
        x1[i] += eps;  // x1: forward perturbation in dimension i
        x2[i] -= eps; // x2: backward perturbation in dimension i
        grad[i] = (evaluate(x1) - evaluate(x2)) / (2 * eps); // Central difference formula
    }
    return grad;
}




std::vector<double> BFGSOptimizer::optimize(const std::vector<double>& x0, int maxIterations) {
    Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(x0.data(), x0.size());
    //hessian matrix  aproximation
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(dim, dim);

    for (int iter = 0; iter < maxIterations; ++iter) {
        Eigen::VectorXd grad = numericalGradient(x);
        if (grad.norm() < 1e-6) break;

        Eigen::VectorXd p = -H * grad; //direction for  (discesa)
        double step = 1e-1; // step a = 0.1 (basic step)

        Eigen::VectorXd x_new = x + step * p; // move to the p direction
        Eigen::VectorXd grad_new = numericalGradient(x_new);

        Eigen::VectorXd s = x_new - x;
        Eigen::VectorXd y = grad_new - grad;

        double sy = s.dot(y);

        //for avoid huge division and calculate new H
        if (sy > 1e-10) {
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim, dim);
            H = (I - s * y.transpose() / sy) * H * (I - y * s.transpose() / sy) + (s * s.transpose()) / sy;
        }

        x = x_new;
    }

    std::vector<double> result(x.data(), x.data() + x.size());
    return result;
}

#endif //USE_EIGEN