#ifndef EUCLIDEANDISTANCE_HPP
#define EUCLIDEANDISTANCE_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class EuclideanDistance : public ObjectiveFunction {
   public:
	virtual double operator()(const std::vector<double>& position) const override;
};

#endif	// EUCLIDEANDISTANCE_HPP