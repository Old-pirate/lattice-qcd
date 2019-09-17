#include "field.h"

Field::Field(unsigned const& length_edge) : Random(length_edge) {}
Field::~Field() {
	this->field_array.clear();
}

void Field::make_field(double const& a) {
	for (size_t i = 0; i < this->size; ++i)
		this->field_array.push_back(a);
}

double Field::calculate_squared_mean_field() {
	double var = 0;
	for (size_t i = 0; i < this->size; ++i)
		var += this->field_array[i];
	return pow(var / this->field_array.size(), 2);
}

void Field::start_time_() {
	this->start_time = clock();
}
void Field::end_time_() {
	this->end_time = clock();
	this->search_time = 0.001*(this->end_time - this->start_time);
}