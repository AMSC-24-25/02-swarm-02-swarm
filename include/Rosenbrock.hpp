#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class Rosenbrock : public ObjectiveFunction {
   public:
	const double a = 1.0;
	const double b = 10.0;

	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// ROSENBROCK_HPP