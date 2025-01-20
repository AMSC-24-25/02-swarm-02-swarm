#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <vector>
#include <cassert>

#include "ObjectiveFunction.hpp"

class Rosenbrock : public ObjectiveFunction {
   public:
	const double a;
	const double b;

	Rosenbrock(const double _a = 1.0, const double _b = 10.0) : a(_a), b(_b) {
		assert(_a >= 0.0);
		assert(_b >= 0.0);
	}

	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// ROSENBROCK_HPP