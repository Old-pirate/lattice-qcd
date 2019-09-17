#include "Random.h"
#include "Field.h"
#include "Ising_model.h"
#include "Scalar_field.h"

int main() {
	setlocale(LC_ALL, "Russian");
	
	{
		Ising_model Im(20);
		ofstream out1("1. модель »зинга - термализаци€ 1.txt", ios_base::trunc);
		ofstream out2("1. модель »зинга - термализаци€ 2.txt", ios_base::trunc);
		Im.definition_of_n_therm(1.0, 500000, out1, out2);
		Im.make_field();
		ofstream out0("1. модель »зинга - фазовый переход.txt", ios_base::trunc);
		Im.phase_transition(1.0, 200000, 2500, out0);
	}
	
	{
		//Scalar_field Sf(5);
		//ofstream out1("2. скал€рное поле - термализаци€ 1.txt", ios_base::trunc);
		//ofstream out2("2. скал€рное поле - термализаци€ 2.txt", ios_base::trunc);
		//Sf.definition_of_n_therm(0.5, 0.01, 0.5, 500000, out1, out2);
		//ofstream out0("2. скал€рное поле - фазовый переход.txt", ios_base::trunc);
		//Sf.make_field();
		//Sf.phase_transition(1.0, 200000, 2500, out0);
	}
	
	cout << endl; system("pause");
	return 0;
}