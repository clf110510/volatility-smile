#ifndef _BS_ANALYTICS
#define _BS_ANALYTICS

#include <cmath>
#include <boost/math/distributions/normal.hpp>

using namespace boost::math;

enum OptionType {Call, Put};
double cnorm(double x);
double bsPricer(OptionType optType, double K, double T, double S_0, double sigma, double r, double q);
double margrabe(double S1, double S2, double T, double q1, double q2, double sigma1, double sigma2, double rho);

double bsPricer(OptionType optType, double K, double T, double F, double sigma);
double bsImpliedVol(OptionType optType, double K, double T, double F, double prc);
   
#endif
