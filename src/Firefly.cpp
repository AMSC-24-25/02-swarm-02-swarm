#include "../include/Firefly.h"

#include <random>

// Constructor: initialize firefly with random position
Firefly::Firefly(int dim) : dimensions(dim), brightness(0.0) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> dis(-50.0, 50.0);	// random range [-5, 5]

	position.resize(dim,0.0);
	// for (int d = 0; d < dim; ++d) {
	// 	position[d] = dis(gen);
	// }
}

// Return current position
std::vector<double> Firefly::getPosition() const {
	return position;
}

// Update position
void Firefly::setPosition(const std::vector<double>& newPosition) {
	position = newPosition;
}

// Return brightness
double Firefly::getBrightness() const {
	return brightness;
}

// Update brightness
void Firefly::setBrightness(double newBrightness) {
	brightness = newBrightness;
}
