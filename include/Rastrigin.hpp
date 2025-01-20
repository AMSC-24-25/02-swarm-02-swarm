#ifndef RASTRIGIN_HPP
#define RASTRIGIN_HPP

#include <vector>
#include <cassert>

#include "ObjectiveFunction.hpp"

class Rastrigin : public ObjectiveFunction {
   public:
	const double A;

	Rastrigin(const double _A = 10.0) : A(_A) {
		assert(_A >= 0.0);
	}

	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// RASTRIGIN_HPP