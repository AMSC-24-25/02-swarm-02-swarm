#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <vector>

#include "ObjectiveFunction.hpp"

class Sphere : public ObjectiveFunction {
   public:
	virtual double operator()(const std::vector<double>& position) const;
};

#endif	// SPHERE_HPP