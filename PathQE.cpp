#include <math.h>
#include "PathQE.h"
#include "RandomProcess.h"
#include "X_simulation.h"
#include "AnalyticPricer.h"

//In this file, we implement the simulation of the variance using the QE Scheme.

//Case where phi <= phi_c : 
//  We compute a and b and use the inverse normal cdf
double scheme_1(double phi, double m, double u)
{
    double b2 = 2 / phi - 1 + sqrt(2 / phi) * sqrt(2 / phi -1);
    double a = m / (1 + b2);
    double Z = inverse_normal_cdf(u);

    return a * pow((Z + sqrt(b2)), 2);

}

//Case where phi > phi_c :
//  We compute p and beta, and we get the volatility by using the standard inverse distribution function method.
double scheme_2(double phi, double m, double u)
{
    double p = (phi - 1)/(phi + 1);
    double beta = (1 - p) / m;
    double inv_distrib = 0;
    if (u > p) {
        inv_distrib = (1 / beta) * log((1 - p)/(1 - u));
    }
    return inv_distrib;
}

//We are computing the next variance having the switching rules and the critical level phi_c
double variance_simulation(double phi_c, double var, double delta, Heston_Model h)
{
    double next_var;
    double m = h.theta + (var - h.theta) * exp(-h.kappa * delta); //m is the mean of the next variance
    double s2 = (1 / h.kappa) * var * h.sigma2 * exp(-h.kappa * delta) * (1 - exp(-h.kappa * delta))
                + h.kappa / 2 * h.theta * h.sigma2 * pow((1 - exp(-h.kappa * delta)), 2); //s2 is the variance of the next variance

    double phi = s2 / pow(m,2);
    double random_number = unif(0, 1);
    if (phi_c < phi) //switching rule
    {
        next_var = scheme_2(phi, m, random_number);
    }
    else
    {
        next_var = scheme_1(phi, m, random_number);
    }
    return next_var;
}

void test_qe_simulator()
{
    class Last : public RandomProcessEval {
    public:
        double virtual operator() (RandomProcess X) const { return X[X.N - 1]; }
    };

    double v0 = .1;
    QE_simulator qe_simulator(10, 1, v0);
    double mean = qe_simulator.monte_carlo_computation(10, Last());
    std::cout << mean << endl;
}

QE_simulator::QE_simulator(int N, double T, double v0, double phi_c, Heston_Model h): x0(v0), phi_c(phi_c), delta(T / N), N(N), T(T), h(h)
{
}

//Simulating a vector of variances having from starting variance
RandomProcess QE_simulator::simulate() const
{   
    double v = x0;
    vector<double> data(N);
    for (int i = 0; i < N; i++) {
        v = variance_simulation(phi_c, v, delta, h);
        data[i] = v;
    }
    return RandomProcess(N, N * delta, data);
}

//Comparing the prices we get using the Broadie Kaya and QE scheme to the closed formula.
void test_convergence_QE(double T, int N, double v0, double x0, double phi_c, double n_mc_simulations){

	Heston_Model hm;
	AnalyticPricer Pricer(T, N, hm);

	QE_simulator* var_simulator_QE = new QE_simulator(N, T, v0, phi_c, hm);
	X_simulator asset_price_simulator_QE(N, T, x0, var_simulator_QE, hm);
	double swap_price = asset_price_simulator_QE.monte_carlo_computation(n_mc_simulations);

    cout << "The price using the closed formula is : " << closed_formula(T, hm, v0) << endl;
    cout << "The price using Broadie Kaya and QE simulation is : " << swap_price << endl;
}
