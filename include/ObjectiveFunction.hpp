#ifndef OBJECTIVEFUNCTION_HPP
#define OBJECTIVEFUNCTION_HPP

#include <vector>

class ObjectiveFunction {
   public:
	virtual double operator()(const std::vector<double>& position) const = 0;

	virtual ~ObjectiveFunction() = default;
};

#endif
