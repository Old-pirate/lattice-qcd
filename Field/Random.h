#pragma once
#include <random>
using namespace std;

class Random
{
protected:
	unsigned length_edge;
	random_device rd;
	mt19937 generator;
	normal_distribution<double>			normal_double;
	uniform_int_distribution<int>		uniform_int;
	uniform_int_distribution<int>		element;
	uniform_real_distribution<double>	random_number_r1;	
	uniform_real_distribution<double>	random_number_r2;

public:
	Random(unsigned const& length_edge);
	virtual ~Random();

	inline int uniform_int_() {
		return this->uniform_int(this->generator) == 0 ? -1 : 1;
	}
};
