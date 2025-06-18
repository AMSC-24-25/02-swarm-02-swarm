#include "State.hpp"


State::State(const std::vector<double>& initialValues, const ObjectiveFunction& objFunc)
    : values(initialValues), objective(objFunc)
{
    cost = objective(initialValues);
}

void State::update(const std::vector<double>& newValues) {
    values = newValues;
    cost = objective(values);
    
    std::cout << std::endl;
}
