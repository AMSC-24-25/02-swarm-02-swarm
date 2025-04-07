#include <vector>
#include <cassert>
#include <limits>

#include "Creature.hpp"

Creature::Creature(const std::vector<double>& position_)
	: position(position_), fitness(std::numeric_limits<double>::infinity()) {
	assert(position_.size() > 0);
}
