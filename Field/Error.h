#pragma once
#include<vector>
#include<iostream>
using namespace std;

class Error
{
	unsigned error_parameter_bin_size;
	double sum_elements, mean_value, sigma;
	vector<double>
		observed_magnitude_series, binned_series, estimated_parameters_series;

public:
	Error();
	~Error();
	Error(Error const&) = delete;
	Error& operator=(Error const&) = delete;

	void make_observed_magnitude_series(double const& var);
	double calculate_sigma(unsigned const& error_parameter_bin_size);

private:
	void make_binned_series();
	void calculate_sum_elements_binned_series();
	void make_estimated_parameters_series();
	void calculate_mean_value_estimated_parameters_series();
	void calculate_sigma_binning_and_jackknife();
};

