#ifndef STATE_HPP
#define STATE_HPP

#include <vector>
#include "ObjectiveFunction.hpp"

class State {
public:
    std::vector<double> values;  // Soluzione
    double cost;                 // Valutazione della funzione obiettivo

    State(const std::vector<double>& initialValues, const ObjectiveFunction& objFunc);

    void update(const std::vector<double>& newValues); // Aggiorna valori + ricostruisce cost

private:
    const ObjectiveFunction& objective;
};

#endif // STATE_HPP
