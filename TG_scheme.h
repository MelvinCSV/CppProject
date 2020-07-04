#pragma once
#include "Simulator.h"
#include "Newton_method.h"
#include "Heston_Model.h"


// Functions for the truncated gaussian scheme 
double fun_m(double var, double delta, Heston_Model h);
double fun_s2(double var, double delta, Heston_Model h);
double fun_r(double psi);
double tg_next_var(double var, double delta, Heston_Model h);

// Simulator for the truncated gaussian scheme
class TG_simulator: public Simulator{};
