#include "Scalar_field.h"

Scalar_field::Scalar_field(unsigned const& length_edge) 
	: Field(length_edge) {
	this->metric = 4;
	this->size = (size_t)pow(this->length_edge, this->metric);
}

void Scalar_field::phase_transition(
	double const& kappa_begin, double const& kappa_end, double const& kappa_step,
	double const& lambda_begin, double const& lambda_end, double const& lambda_step, 
	double const& epsilon, unsigned const& n_therm,
	unsigned const& error_parameter_bin_size, ostream& out) {
	unsigned x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0;
	
	cout << "Phase transition\n";
	out << "lambda\t" << "kappa\t" << "Squared mean scalar field\t" << "sigma"
		<< endl;

	for (double lambda = lambda_begin; lambda <= lambda_end; lambda += lambda_step) {
		start_time_();
		for (double kappa = kappa_begin; kappa < kappa_end; kappa += kappa_step) {
			thermalization(kappa, lambda, epsilon, n_therm, x_1, x_2, x_3, x_4);
			out << left
				<< setw(6) << lambda << "\t"
				<< setw(6) << kappa << "\t"
				<< setw(25) << calculate_squared_mean_field() << "\t"
				<< setw(10) << err.calculate_sigma(error_parameter_bin_size)
				<< endl;
		}
		out << endl;

		end_time_(); 
		cout
			<< "\tЗатраченное время на вычисления kappa[0.01, 0.50] при lambda = " 
			<< lambda << " --- " << this->search_time << " с" << endl;
	}
}
void Scalar_field::thermalization(
	double const& kappa
	, double const& lambda
	, double const& epsilon
	, unsigned const& n_therm
	, unsigned& x_1
	, unsigned& x_2
	, unsigned& x_3
	, unsigned& x_4
)
{
	for (int i = 0; i <= n_therm; ++i) {
		choose_node_metropolis_algorithm(kappa, lambda, epsilon, x_1, x_2, x_3, x_4);
		err.make_observed_magnitude_series(calculate_squared_mean_field());
	}
}

void Scalar_field::definition_of_n_therm(
	double const& kappa, 
	double const& lambda, 
	double const& epsilon,
	unsigned const& n_therm,
	ostream& out1, 
	ostream& out2
) {
	make_field(1); 
	start_time_();
	thermalization(kappa, lambda, epsilon, n_therm, out1);
	end_time_();
	cout << "Затраченное время на выполнение definition_of_n_therm.part.1 " 
		 << this->search_time << " с" << endl;

	this->field_array.clear();
	make_field(0); 
	start_time_();
	thermalization(kappa, lambda, epsilon, n_therm, out2);
	end_time_();
	cout << "Затраченное время на выполнение definition_of_n_therm.part.2 "
		 << this->search_time << " с" << endl;
	
	this->field_array.clear();
}
void Scalar_field::thermalization(
	double const& kappa, double const& lambda, double const& epsilon,
	unsigned const& n_therm, ostream& out) {
	unsigned x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0;
	for (int i = 0; i <= n_therm; ++i) {
		out << i << "\t"
			<< calculate_squared_mean_field() << "\t"
			<< calculate_action(kappa, lambda)
			<< endl;
		choose_node_metropolis_algorithm(kappa, lambda, epsilon, x_1, x_2, x_3, x_4);
	}
}

void Scalar_field::definition_of_epsilon(
	double const& kappa, double const& lambda, 
	unsigned const& n_therm, ostream& out) {
	start_time_();
	for (int i = 0; i <= 5; ++i) {
		for (double epsilon = 0; epsilon <= 2.1; epsilon += 0.1) {
			out << epsilon << "\t" 
				<< calculate_overlap_factor(kappa, lambda, epsilon, n_therm) 
				<< endl;
		}
		out << endl;
	}
	end_time_();
	cout << "Затраченное время на выполнение definition_of_epsilon"
		 << this->search_time << " с" << endl;
}
double	Scalar_field::calculate_overlap_factor(
	double const& kappa, double const& lambda, double const& epsilon, 
	unsigned const& n_therm) {
	unsigned  x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0;
	double var = 0;

	Field::make_field(1);
	for (int i = 0; i <= n_therm; ++i) {
		x_1 = this->element(this->generator); x_2 = this->element(this->generator);
		x_3 = this->element(this->generator); x_4 = this->element(this->generator);
	
		var += metropolis_algorithm(kappa, lambda, epsilon, x_1, x_2, x_3, x_4);
	}
	this->field_array.clear();

	return var / n_therm * 100;
}

void Scalar_field::choose_node_metropolis_algorithm(
	double const& kappa, double const& lambda, double const& epsilon,
	unsigned& x_1, unsigned& x_2, unsigned& x_3, unsigned& x_4) {
	x_1 = this->element(this->generator); x_2 = this->element(this->generator);
	x_3 = this->element(this->generator); x_4 = this->element(this->generator);

	metropolis_algorithm(kappa, lambda, epsilon, x_1, x_2, x_3, x_4);
}
unsigned Scalar_field::metropolis_algorithm(
	double const& kappa
	, double const& lambda
	, double const& epsilon
	,	unsigned const& x_1
	, unsigned const& x_2
	, 	unsigned const& x_3
	, unsigned const& x_4
) {
	double
		old_spin = this->field_array[
			pow(this->length_edge, 0)*(x_4)+
			pow(this->length_edge, 1)*(x_3)+
			pow(this->length_edge, 2)*(x_2)+
			pow(this->length_edge, 3)*(x_1)],
		new_spin = old_spin + epsilon*this->random_number_r2(this->generator);
	
	field_array[
		pow(this->length_edge, 0)*(x_4) +
		pow(this->length_edge, 1)*(x_3) +
		pow(this->length_edge, 2)*(x_2) +
		pow(this->length_edge, 3)*(x_1)] = old_spin;
	double old_Action = calculate_action(kappa, lambda);

	field_array[
		pow(this->length_edge, 0)*(x_4) +
		pow(this->length_edge, 1)*(x_3) +
		pow(this->length_edge, 2)*(x_2) +
		pow(this->length_edge, 3)*(x_1)] = new_spin;
	double new_Action = calculate_action(kappa, lambda);

	double delta = new_Action - old_Action;

	if (delta >= 0) {
		double r = this->random_number_r1(this->generator);

		if (r <= exp(-delta)) {
			field_array[
				pow(this->length_edge, 0)*(x_4) +
				pow(this->length_edge, 1)*(x_3) +
				pow(this->length_edge, 2)*(x_2) +
				pow(this->length_edge, 3)*(x_1)] = new_spin;
			return 1;
		}
		else {
			field_array[
				pow(this->length_edge, 0)*(x_4) +
				pow(this->length_edge, 1)*(x_3) +
				pow(this->length_edge, 2)*(x_2) +
				pow(this->length_edge, 3)*(x_1)] = old_spin;
			return 0;
		}
	}
	else {
		field_array[
			pow(this->length_edge, 0)*(x_4) +
			pow(this->length_edge, 1)*(x_3) +
			pow(this->length_edge, 2)*(x_2) +
			pow(this->length_edge, 3)*(x_1)] = new_spin;
		return 1;
	}
	return 0;
}
double Scalar_field::calculate_action(double const& kappa, double const& lambda) {
	double A, B, C, var = 0;
	for (int x_4 = 0; x_4 < this->length_edge; ++x_4) {
		for (int x_3 = 0; x_3 < this->length_edge; ++x_3) {
			for (int x_2 = 0; x_2 < this->length_edge; ++x_2) {
				for (int x_1 = 0; x_1 < this->length_edge; ++x_1) {
					A = this->field_array[
							pow(this->length_edge, 0)*(x_4) +
							pow(this->length_edge, 1)*(x_3) +
							pow(this->length_edge, 2)*(x_2) +
							pow(this->length_edge, 3)*(x_1)] * (
						this->field_array[
							pow(this->length_edge, 0)*(x_4) +
							pow(this->length_edge, 1)*(x_3) +
							pow(this->length_edge, 2)*(x_2) +
							pow(this->length_edge, 3)*((x_1 + 1) % this->length_edge)] +
						this->field_array[
							pow(this->length_edge, 0)*(x_4) +
							pow(this->length_edge, 1)*(x_3) +
							pow(this->length_edge, 2)*((x_2 + 1) % this->length_edge) + 
							pow(this->length_edge, 3)*(x_1)] +
						this->field_array[
							pow(this->length_edge, 0)*(x_4) +
							pow(this->length_edge, 1)*((x_3 + 1) % this->length_edge) + 
							pow(this->length_edge, 2)*(x_2) +
							pow(this->length_edge, 3)*(x_1)] +
						this->field_array[
							pow(this->length_edge, 0)*((x_4 + 1) % this->length_edge) +
							pow(this->length_edge, 1)*(x_3) +
							pow(this->length_edge, 2)*(x_2) +
							pow(this->length_edge, 3)*(x_1)]);
							B = pow(this->field_array[pow(this->length_edge, 0)*(x_4)
								+pow(this->length_edge, 1)*(x_3)
								+pow(this->length_edge, 2)*(x_2)
								+pow(this->length_edge, 3)*(x_1)
							], 2);
					C = pow(
						pow(this->field_array[
							pow(this->length_edge, 0)*(x_4) +
							pow(this->length_edge, 1)*(x_3) +
							pow(this->length_edge, 2)*(x_2) +
							pow(this->length_edge, 3)*(x_1)], 2)
						- 1, 2);
					var += -2.0 * kappa*A + B + lambda*C;
				}
			}
		}
	}
	return (var - lambda);
}

void Scalar_field::make_field() {
	for (size_t i = 0; i < this->size; ++i)
		this->field_array.push_back(this->normal_double(this->generator));
}