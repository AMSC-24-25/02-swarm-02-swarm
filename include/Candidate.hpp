#ifndef CANDIDATE_HPP
#define CANDIDATE_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

#include "ObjectiveFunction.hpp"

class Candidate {
	public:
	double f0;
	std::vector<double> candidate;

	private:
	const ObjectiveFunction& func;

	public:
	Candidate(size_t dimensions_, double lower_bound_, double upper_bound_, size_t seed_, const ObjectiveFunction& func_);
	Candidate(const std::vector<double>& candidate_, const ObjectiveFunction& func);
	void updatePosition(const std::vector<double>& position);

};

#endif