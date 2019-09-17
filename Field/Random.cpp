#include"Random.h"
#include<iostream>

Random::Random(unsigned const& length_edge)
	: length_edge(length_edge), generator(rd())
	, normal_double(0, 1), uniform_int(0, 1), element(0, length_edge - 1)
	, random_number_r1(0, 1), random_number_r2(-1, 1) {
}
Random::~Random() {}