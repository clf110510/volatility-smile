#ifndef _IMPLIEDVOL_H
#define _IMPLIEDVOL_H

#include <vector>
#include <utility>
#include <memory>
#include "BSAnalytics.h"

using namespace std;

class Smile
{
 public:
  // constructor
  Smile( const vector< pair<double, double> >& _marks ) : marks(_marks) {}
  
  virtual void Construct() = 0;
  virtual double Vol(double strike) const = 0;

 protected:
  // strike to implied vol marks
  vector< pair<double, double> > marks;
  
};

// CubicSpline interpolated smile, extrapolate flat
class CubicSplineSmile : public Smile
{
 public:
  // constructor
  CubicSplineSmile( const vector< pair<double, double> >& _marks ) : Smile(_marks) {}
  
  virtual void Construct();
  virtual double Vol(double strike) const;

 private:
  vector<double> y2; // second derivatives
  
};

// Arbitrage free interpolated smile, extrapolate flat
class ArbitrageFreeSmile : public Smile
{
 public:
  ArbitrageFreeSmile(const vector<pair<double,double>>& _marks) : Smile(_marks) {}
  ArbitrageFreeSmile(const vector<pair<double,double>>& _marks,
               double _F, double _T)
    : Smile(_marks), F(_F), T(_T) {}

  virtual void Construct();
  virtual double Vol(double strike) const;

 protected:
  int SearchStrike(double strike) const;

 private:
  vector<double> K;
  vector<double> C;
  vector<double> y2;
  double F, T;
};

class ImpliedVol 
{
 public:
  ImpliedVol(const vector<pair<double, shared_ptr<Smile>>>&);

  // linear interpolation in variance, along the strike line
  double Vol(double t, double k) const;

  double dVoldK(double t, double k) const {
    return (Vol(t, k+0.00001) - Vol(t, k-0.00001)) / 0.00002;
  }

  double dVoldT(double t, double k) const {
    return (Vol(t+0.0005, k) - Vol(t, k)) / 0.0005;
  }

  double dVol2dK2(double t, double k) const {
    return (Vol(t, k+0.00001) + Vol(t, k-0.00001) -  2*Vol(t, k))
           / 0.00001 / 0.00001;
  }
  
 private:
  vector<pair<double, shared_ptr<Smile>>> pillarSmiles;
};

#endif
