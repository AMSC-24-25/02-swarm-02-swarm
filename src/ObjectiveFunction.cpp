#include "ObjectiveFunction.hpp"


ObjectiveFunction::ObjectiveFunction(){
    sum = 0.0;
}

double ObjectiveFunction::getValueFunction(const vector<double>& position){
    sum = 0.0;
    for (double x : position) {
        sum += x * x;
    }
    return sum;
}