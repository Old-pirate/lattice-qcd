#pragma once
#include "Random.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <ctime>

class Field : public Random
{
protected:
	vector<double> field_array;
	size_t size;
	unsigned metric;
	double start_time, end_time, search_time;

public:
	Field(unsigned const& length_edge);
	virtual ~Field();
	Field(Field const&) = delete;
	Field & operator=(Field const&) = delete;
	
	friend ostream& operator<<(ostream& out, const Field& dt) {
		out 
			<< "length_edge = " << dt.length_edge << "; metric = " 
			<< dt.metric << "; size = " << dt.size << ";";
		return out;
	}
	double & operator[](size_t i) {
		return this->field_array[i];
	}
	
	virtual void make_field() = 0;
	
	void make_field(double const& a);
	
	double calculate_squared_mean_field();

	void start_time_();
	void end_time_();
};