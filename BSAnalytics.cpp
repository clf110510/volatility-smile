#include "BSAnalytics.h"

double cnorm(double x)
{
  return cdf(normal(), x);
  
  // constants --- implementation without boost and erf
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  return 0.5*(1.0 + sign*y);
}

double bsPricer(OptionType optType, double K, double T, double S_0, double sigma, double r, double q)
{
  double sigmaSqrtT = sigma * std::sqrt(T);
  double d1 = (std::log(S_0 / K) + (r - q)*T)/sigmaSqrtT + 0.5 * sigmaSqrtT;
  double d2 = d1 - sigmaSqrtT;

  double V_0;
  switch (optType)
    {
    case Call:
      V_0 = S_0 * exp(-q*T) * cnorm(d1) - K * exp(-r*T) * cnorm(d2);
      break;
    case Put:
      V_0 = K * exp(-r*T) * cnorm(-d2) - S_0 * exp(-q*T) * cnorm(-d1);
      break;
    default:
      throw "unsupported optionType";
    }
  return V_0;
}

double bsPricer(OptionType optType, double K, double T, double F, double sigma)
{
  double sigmaSqrtT = sigma * std::sqrt(T);
  double d1 = std::log(F / K) / sigmaSqrtT + 0.5 * sigmaSqrtT;
  double d2 = d1 - sigmaSqrtT;

  double V_0;
  switch (optType)
    {
    case Call:
      V_0 = F * cnorm(d1) - K * cnorm(d2);
      break;
    case Put:
      V_0 = K * cnorm(-d2) - F * cnorm(-d1);
      break;
    default:
      throw "unsupported optionType";
    }
  return V_0;
}

double bsImpliedVol(OptionType optType, double K, double T, double F, double prc)
{
  const double eps = 1e-6;

  // binary search
  double a = eps;
  double b = 10 + eps;

  // in case of the left side is not small enough
  if (bsPricer(optType, K, T, F, a) > prc) {
    while (1) {
      b = a;
      a /= 2;
      if (bsPricer(optType, K, T, F, a) < prc) {
        break;
      }
    }
  }

  // in case of the right side is not large enough
  if (bsPricer(optType, K, T, F, b) < prc) {
    while(1) {
      a = b;
      b *= 2;
      if (bsPricer(optType, K, T, F, b) > prc) {
        break;
      }
    }
  }

  // start binary search
  while (1) {
    double c = (a + b) / 2;
    double new_prc = bsPricer(optType, K, T, F, c);
    if (std::abs(new_prc - prc) / prc < eps) {
      return c;
    } else {
      if (new_prc > prc) {
        b = c;
      } else {
        a = c;
      }
    }
  }
}

double margrabe(double S1, double S2, double T, double q1, double q2, double sigma1, double sigma2, double rho)
{
  double sigma = sqrt(pow(sigma1, 2) + pow(sigma2, 2) - 2 * sigma1*sigma2*rho);
  double d1 = (log(S1 / S2) + (q2 - q1 + pow(sigma, 2) / 2)*T) / (sigma*sqrt(T));
  double d2 = d1 - sigma*sqrt(T);  
  return exp(-q1*T)*S1*cnorm(d1) - exp(-q2*T)*S2*cnorm(d2);
}
