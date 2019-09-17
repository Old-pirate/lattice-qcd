#pragma once
#include "Field.h"
#include "Error.h"

class Scalar_field :
	public Field
{
private:
	Error err;
public:
	using Field::make_field;

	Scalar_field(unsigned const& length_edge);
	Scalar_field(Scalar_field const&) = delete;
	Scalar_field & operator=(Scalar_field const&) = delete;

	void phase_transition(
		double const& kappa_begin, double const& kappa_end, double const& kappa_step,
		double const& lambda_begin, double const& lambda_end, double const& lambda_step,
		double const& epsilon, unsigned const& n_therm,
		unsigned const& error_parameter_bin_size, ostream& out);
	void thermalization(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned const& n_therm,
		unsigned& x_1, unsigned& x_2, unsigned& x_3, unsigned& x_4);

	void definition_of_n_therm(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned const& n_therm, ostream& out1, ostream& out2);
	void thermalization(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned const& n_therm, ostream& out);
	
	void definition_of_epsilon(
		double const& kappa, double const& lambda,
		unsigned const& n_therm, ostream& out);
	double calculate_overlap_factor(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned const& n_therm);

	void choose_node_metropolis_algorithm(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned& x_1, unsigned& x_2, unsigned& x_3, unsigned& x_4);
	unsigned metropolis_algorithm(
		double const& kappa, double const& lambda, double const& epsilon,
		unsigned const& x_1, unsigned const& x_2, 
		unsigned const& x_3, unsigned const& x_4);
	double calculate_action(double const& kappa, double const& lambda);
	
	void make_field();
};

