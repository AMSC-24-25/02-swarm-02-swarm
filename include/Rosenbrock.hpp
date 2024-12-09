#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class Rosenbrock : public ObjectiveFunction {
   public:
	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// ROSENBROCK_HPP