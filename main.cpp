#include "AnalyticPricer.h"
#include <iostream>
#include "RandomProcess.h"
#include <iostream>
#include "PathQE.h"
#include "X_simulation.h"
#include "Newton_method.h"
#include <vector>
#include<time.h>


using namespace std;

int main()
{	
	srand(time(NULL));

	// Randomize the seed

	//double T = 10;
	//int N = 1000;
	//Heston_Model hm;
	//AnalyticPricer Pricer(T, N, hm);
	//double v0 = theta0 * 1.5;
	//double x0 = 1;
	//double phi_c = 1.5;
	//int n_mc_simulations = 10;

	//cout << Pricer.compute_swap_price(v0) << endl;
	//cout << closed_formula(T, hm, v0) << endl;

	//QE_simulator* var_simulator_QE = new QE_simulator(N, T, v0, phi_c, hm);
	//X_simulator asset_price_simulator_QE(N, T, x0, var_simulator_QE, hm);
	//double swap_price = asset_price_simulator_QE.monte_carlo_computation(n_mc_simulations);

	//cout << swap_price << endl;

	test_qe_simulator();
	//double x_test = X_next_simulate(1, 0, 0, 1, 0.05, 0, 0.5, 1);
	//double var_test = variance_simulation(1.5, 0.03, 0.5, 0.01, 1, 0.1);
	//cout << x_test << endl;
	//cout << var_test << endl;

	//double y = test_newton_method();
	//cout << y;
	
	//int N = 10;
	//vector<double> v(N);
	//for (int i = 0; i < N; i++)
	//	v[i] = 1;

	//cout << v[N - 1] << endl;

	return 0;
}

