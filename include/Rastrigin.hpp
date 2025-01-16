#ifndef RASTRIGIN_HPP
#define RASTRIGIN_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class Rastrigin : public ObjectiveFunction {
   public:
	const double A = 10.0;

	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// RASTRIGIN_HPP