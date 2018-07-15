#include "LocalVol.h"
#include <cmath>
#include <iostream>

double LocalVol::LV(double t, double s) const
{
  double imp = iv.Vol(t, s);
  if (t < 1e-6)
      return imp;
  double dvdk = iv.dVoldK(t, s);
  double dvdt = iv.dVoldT(t, s);
  double d2vdk2 = iv.dVol2dK2(t, s);
  double d1 = (std::log(S0/s) + (rd-rf)*t + 0.5 * imp * imp * t ) / imp / std::sqrt(t);
  double nominator = imp*imp + 2 * t * imp * dvdt + 2 * (rd - rf) * s * t * imp * dvdk;
  double denominator = std::pow(1+s*d1*sqrt(t)*dvdk, 2) + s * s * t * imp * (d2vdk2 - d1 * sqrt(t) * dvdk * dvdk);
  // double denominator = std::pow(1+s*d1*t*dvdk, 2) + s * s * t * imp * (d2vdk2 - d1 * t * dvdk * dvdk);
  // double d2 = (std::log(S0/s) + (rd-rf)*t - 0.5 * imp * imp * t ) / imp / std::sqrt(t);
  // double denominator = 1 + 2* s * d1 * sqrt(t) * dvdk + s * s * t * (d1 * d2 * dvdk * dvdk + imp * d2vdk2);
  double localvar = std::min(std::max(nominator / denominator, 1e-8), 1.0);  
  if (nominator <= 0)
    localvar = 1e-8;
  if (denominator <= 0)
    localvar = 1.0;
  return std::sqrt(localvar);
}
