#pragma once
#include "Field.h"
#include "Error.h"

class Ising_model : public Field
{
	Error err;
public:
	using Field::make_field;

	Ising_model(unsigned const& length_edge);
	
	void phase_transition(double const& T_begin, unsigned const& n_therm,
		unsigned const& error_parameter_bin_size, ostream& out);
	void thermalization(double const& T, unsigned const& n_therm, 
		unsigned& x_1, unsigned& x_2);

	void definition_of_n_therm(double const& T, unsigned const& n_therm, 
		ostream& out1, ostream& out2);
	void thermalization(double const& T, unsigned const& n_therm, ostream &out);

	void choose_node_metropolis_algorithm(double const& T, 
		unsigned& x_1, unsigned& x_2);
	void metropolis_algorithm(double const& T, unsigned const& x_1, unsigned const& x_2);
	int calculate_hamiltonian(unsigned const& x_1, unsigned const& x_2);
	
	void make_field();
};

