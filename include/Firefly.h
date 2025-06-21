#ifndef FIREFLY_H
#define FIREFLY_H
#include <vector>

class Firefly {
public:
    //  Costruttore
    Firefly(int dimensions);

    //  Ottieni la posizione della lucciola
    std::vector<double> getPosition() const;

    //  Imposta la posizione della lucciola
    void setPosition(const std::vector<double>& newPosition);

    //  Ottieni la luminosità della lucciola
    double getBrightness() const;

    //  Imposta la luminosità della lucciola
    void setBrightness(double newBrightness);

private:
    int dimensions;  //  Dimensione dello spazio di ricerca
    std::vector<double> position;  //  Posizione della lucciola
    double brightness;  //  Luminosità della lucciola
};




#endif //FIREFLY_H
