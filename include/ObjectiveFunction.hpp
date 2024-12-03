#ifndef OBJECTIVEFUNCTION_HPP
#define OBJECTIVEFUNCTION_HPP

#include <vector>

class ObjectiveFunction {
   public:
	double sum;
	ObjectiveFunction();

	double getValueFunction(const std::vector<double>& position);
};

#endif
