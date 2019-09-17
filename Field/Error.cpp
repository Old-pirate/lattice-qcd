#include "Error.h"

Error::Error() {}
Error::~Error() {
	observed_magnitude_series.clear();
	binned_series.clear();
	estimated_parameters_series.clear();
}

void Error::make_observed_magnitude_series(double const& var) {
	this->observed_magnitude_series.push_back(var);
}

double Error::calculate_sigma(unsigned const& error_parameter_bin_size) {
	this->error_parameter_bin_size = error_parameter_bin_size;
	
	make_binned_series();
	observed_magnitude_series.clear();
	
	calculate_sum_elements_binned_series();
	make_estimated_parameters_series();
	binned_series.clear();

	calculate_mean_value_estimated_parameters_series();
	calculate_sigma_binning_and_jackknife();
	estimated_parameters_series.clear();

	return this->sigma;
}

void Error::make_binned_series() {
	double var = 0;
	for (size_t i = 0; i < this->observed_magnitude_series.size(); ++i) {
		var += this->observed_magnitude_series[i];
		if ((i + 1) % this->error_parameter_bin_size == 0) {
			this->binned_series.push_back(var / this->error_parameter_bin_size);
			var = 0;
		}
	}
}
void Error::calculate_sum_elements_binned_series() {
	double var = 0;
	for (size_t i = 0; i < this->binned_series.size(); ++i)
		var += this->binned_series[i];
	this->sum_elements = var;
}
void Error::make_estimated_parameters_series() {
	double estimated_parameter = 0; 
	for (size_t i = 0; i < this->binned_series.size(); ++i) {
		estimated_parameter = (this->sum_elements - this->binned_series[i]) /
			(this->binned_series.size() - 1);
		this->estimated_parameters_series.push_back(estimated_parameter);
	}
}
void Error::calculate_mean_value_estimated_parameters_series() {
	double var = 0;
	for (size_t i = 0; i < this->estimated_parameters_series.size(); ++i)
		var += this->estimated_parameters_series[i];
	this->mean_value = var / this->estimated_parameters_series.size();
}
void Error::calculate_sigma_binning_and_jackknife() {
	double var = 0;
	for (size_t i = 0; i < this->estimated_parameters_series.size(); ++i)
		var += pow(this->estimated_parameters_series[i] - this->mean_value, 2);
	this->sigma = sqrt(var * (this->estimated_parameters_series.size() - 1) /
		this->estimated_parameters_series.size());
}
