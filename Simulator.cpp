#include "Simulator.h"

double Simulator::monte_carlo_computation(int n_simulation, RandomProcessEval const& f)
{	
	double res = 0;
	RandomProcess X;
	for (int i = 0; i < n_simulation; i++) {
		X = this->simulate();
		res += f(X);
	}

	return res / n_simulation;
}
