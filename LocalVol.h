#ifndef _LOCALVOL_H
#define _LOCALVOL_H

#include "ImpliedVol.h"

using namespace std;

class LocalVol
{
 public:
  LocalVol( const ImpliedVol & _iv, double _S0, double _rd, double _rf)
   : iv(_iv), rd(_rd), rf(_rf), S0(_S0) {}

  double LV(double t, double s) const;

 private:
  ImpliedVol iv;
  double S0;
  double rd;
  double rf;
};

#endif 
