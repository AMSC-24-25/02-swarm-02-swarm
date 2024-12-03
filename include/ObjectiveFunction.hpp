#ifndef OBJECTIVEFUNCTION_HPP
#define OBJECTIVEFUNCTION_HPP

#include <vector>

using namespace std;

class ObjectiveFunction{
    public:
    double sum;
    ObjectiveFunction();

    double getValueFunction(const vector<double>& position);

};

#endif