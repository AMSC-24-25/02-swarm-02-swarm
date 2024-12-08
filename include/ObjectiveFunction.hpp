#ifndef OBJECTIVEFUNCTION_HPP
#define OBJECTIVEFUNCTION_HPP

#include <vector>

class ObjectiveFunction {
   public:
	double operator()(const std::vector<double>& position) const;
};

#endif
