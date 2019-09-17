#include "ising_model.h"

Ising_model::Ising_model(unsigned const& length_edge) 
	: Field(length_edge) {
	this->metric = 2;
	this->size = (size_t)pow(this->length_edge, this->metric);
}


void Ising_model::phase_transition(double const& T_begin, unsigned const& n_therm,
	unsigned const& error_parameter_bin_size, ostream& out) {
	unsigned x_1 = 0, x_2 = 0;

	start_time_();
	for (double T = T_begin; T <= 3.1; T += 0.1) {
		thermalization(T, n_therm, x_1, x_2);
		out << left 
			<< setw(10) << T 
			<< setw(16) << calculate_squared_mean_field() 
			<< setw(16) << err.calculate_sigma(error_parameter_bin_size)
			<< endl;
	}
	end_time_();
	cout << "Затраченное время на выполнение phase_transition " 
		<< this->search_time << " с" << endl;
}
void Ising_model::thermalization(double const& T, unsigned const& n_therm,
	unsigned& x_1, unsigned& x_2) {
	for (int i = 0; i <= n_therm; ++i) {
		choose_node_metropolis_algorithm(T, x_1, x_2);
		err.make_observed_magnitude_series(calculate_squared_mean_field());
	}
}

void Ising_model::definition_of_n_therm(double const& T, unsigned const& n_therm, 
	ostream& out1, ostream& out2) {
	make_field(1);
	start_time_();
	thermalization(T, n_therm, out1);
	end_time_();
	cout << "Затраченное время на выполнение definition_of_n_therm.part.1 " 
		 << this->search_time << " с" << endl;

	this->field_array.clear();
	make_field(0); 
	start_time_();
	thermalization(T, n_therm, out2);
	end_time_(); 
	cout << "Затраченное время на выполнение definition_of_n_therm.part.2 " 
		 << this->search_time << " с" << endl;
}
void Ising_model::thermalization(double const &T, unsigned const& n_therm,
	ostream& out) {
	unsigned x_1 = 0, x_2 = 0;
	for (int i = 0; i <= n_therm; ++i) {
		out << calculate_squared_mean_field() << "\t" << i << endl;
		choose_node_metropolis_algorithm(T, x_1, x_2);
	}
}

void Ising_model::choose_node_metropolis_algorithm(double const& T,
	unsigned& x_1, unsigned& x_2) {
	x_1 = this->element(this->generator);
	x_2 = this->element(this->generator);

	metropolis_algorithm(T, x_1, x_2);
}
void Ising_model::metropolis_algorithm(double const& T,
	unsigned const& x_1, unsigned const& x_2) {
	int	old_spin = this->field_array[x_2 + this->length_edge*x_1];
	int old_hamiltonian = calculate_hamiltonian(x_1, x_2);

	int new_spin = uniform_int_();
	this->field_array[x_2 + this->length_edge*x_1] = new_spin;
	int new_hamiltonian = calculate_hamiltonian(x_1, x_2);

	int delta = new_hamiltonian - old_hamiltonian;

	if (delta >= 0) {
		double r = this->random_number_r1(this->generator);

		if (r <= exp(-delta / T))
			this->field_array[x_2 + this->length_edge*x_1] = new_spin;
		else
			this->field_array[x_2 + this->length_edge*x_1] = old_spin;
	}
	else
		this->field_array[x_2 + this->length_edge*x_1] = new_spin;
}
int Ising_model::calculate_hamiltonian(unsigned const& x_1, unsigned const& x_2) {
	int var = this->field_array[x_2 + this->length_edge*x_1] * (
		this->field_array[(x_2 - 1 + this->length_edge) % this->length_edge + this->length_edge*x_1] +
		this->field_array[(x_2 + 1) % this->length_edge + this->length_edge*x_1] +
		this->field_array[x_2 + this->length_edge*((x_1 - 1 + this->length_edge) % this->length_edge)] +
		this->field_array[x_2 + this->length_edge*((x_1 + 1) % this->length_edge)]);
	return -var;
}

void Ising_model::make_field() {
	for (size_t i = 0; i < this->size; ++i)
		this->field_array.push_back(uniform_int_());
}
