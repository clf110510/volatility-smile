#include "ImpliedVol.h"
#include "LocalVol.h"
#include "BSAnalytics.h"
#include "Solver.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <random>

#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric;
using namespace std;

// we take number of paths and number of time steps as input
vector< pair<double,double> >
mcLocalVol( function<vector<double> (double)> payoff, int nPaths, int nT
	  , double S, double T, double r, double q, const LocalVol& lv)
{
  // Generate a normal distribution
  std::seed_seq seed{5};
  std::mt19937 e(seed);
  std::normal_distribution<> normal_dist(0, 1);
  
  std::vector<double> sum;
  std::vector<double> hsquare;
  double dt = T / nT;
  for(int i = 0; i < nPaths; i++) {
    // for each path, we simulate until maturity T
    double X = std::log(S);
    for(int j = 1; j <= nT; j++) {
      double sigma = lv.LV((j-1)*dt, std::exp(X)); // query local vol
      // std::cout << (j-1)*dt << "," << std::exp(X) << " = " << sigma << std::endl;
      double a = (r - q - 0.5 * sigma * sigma) * dt; // drift term
      double b = normal_dist(e) * sqrt(dt) * sigma;  // diffusion term
      X += a + b; // update state variable
    }

    // collect the sum and sum sqaure of h
    std::vector<double> hs = payoff(std::exp(X));
    if ( i == 0 ) {
        sum = hs;
        for ( int k = 0; k < hs.size(); k++ ) {
            hsquare.push_back(hs[k] * hs[k]);
        }
    }
    else {
        for ( int k = 0; k < sum.size(); k++ ) {
            sum[k] += hs[k];
            hsquare[k] += hs[k] * hs[k];
        }
    }
  }
  vector< pair<double, double> > res(sum.size());
  for ( int i = 0; i < sum.size(); i++ )
  {
      res[i].first = std::exp(-r*T) * sum[i] / nPaths;
      res[i].second = sqrt((hsquare[i] / nPaths - (sum[i]/nPaths)*(sum[i]/nPaths))/nPaths);
  }
  return res;
}

// market strangle to smile strangle, note that 
vector< pair<double, double> > input2Marks ( double spot, double rd, double rf, double t, double atmvol,
                                             double ms25, double rr25, double ms10, double rr10 )
{
    // first we need to convert market strangle to smile strangle (butterfly, or bf)
    double bf25 = ms25, bf10 = ms10; // ommitted the conversion, just take the market strangle as smile strangle.

    // this below takes all the input and calculates smile strangle at 25 and 10 delta
    // then resolve the marks in (strike, vol) pairs
    double atmStrike = spot * exp((rd-rf)*t); // ATMForward convention
    double c25 = atmvol + bf25 + rr25 * 0.5;  // call25 implied vol
    double p25 = atmvol + bf25 - rr25 * 0.5;  // put25 implied vol
    double c10 = atmvol + bf10 + rr10 * 0.5;  // call25 implied vol
    double p10 = atmvol + bf10 - rr10 * 0.5;  // put25 implied vol
    
    // we can use the function quantile to avoid calling a root searcher
    // assuming forward delta convention, call delta is N(d1) and put delta is N(d1) - 1
    boost::math::normal norm(0, 1.0); 
    double fwd = atmStrike;
    double c25_strike = fwd / exp((quantile(norm, 0.25) * c25 * sqrt(t) - 0.5 * c25 * c25 * t));
    double c10_strike = fwd / exp((quantile(norm, 0.10) * c10 * sqrt(t) - 0.5 * c10 * c10 * t));
    double p25_strike = fwd / exp((quantile(norm, 0.75) * p25 * sqrt(t) - 0.5 * p25 * p25 * t));
    double p10_strike = fwd / exp((quantile(norm, 0.90) * p10 * sqrt(t) - 0.5 * p10 * p10 * t));
    
    // we may want to check if these strikes are ascending --- omitted again
    vector< pair<double, double> > marks(5);
    marks[0] = pair<double, double>(p10_strike, p10);
    marks[1] = pair<double, double>(p25_strike, p25);
    marks[2] = pair<double, double>(atmStrike, atmvol);
    marks[3] = pair<double, double>(c25_strike, c25);
    marks[4] = pair<double, double>(c10_strike, c10);
    return marks; 
}

double StrikeFromDelta(double spot, double rd, double rf, double delta, double t, const ImpliedVol& iv)
{
    double fwd = spot * exp((rd - rf) * t);
    auto fun = [fwd, delta, iv, t](double k)
               {
                   double sigma = iv.Vol(t, k);
                   double d1 = (log(fwd / k) + 0.5 * sigma * sigma * t) / sigma / sqrt(t);
                   return cnorm(d1) - delta;
               };
    double lowerBound = 0.001, upperBound = 1.0;
    RootBracketing(fun, lowerBound, upperBound, 0.00001);
    return rfbrent(fun, lowerBound, upperBound, 1e-6);
}

// output implied vol and local vol to txt files for illustration or plot using testPlot.py
void ShowIVandLV(const ImpliedVol& iv, const LocalVol & lv, double kMin, double kMax)
{
  ofstream fout3("impliedvol.txt");
  ofstream fout4("localvol.txt");
  double dk = (kMax - kMin) / 50;
  double dt = 10 / 60.0;
  for (int i = 0; i < 60; i++) {
    double t = 0.02 + dt * i;
    for (int j = 0; j < 50; j++) {
      double k = kMin + dk * j;
      fout3 << t << "\t" << k << "\t " << iv.Vol(t, k) << endl;
      fout4 << t << "\t" << k << "\t " << lv.LV(t, k) << endl;
    }
  }
}

// 2D array is very difficult to use --- double v[][19], better to use matrix<double>
void ShowArray2(double v[][19], double tenors[], double nDeltas, double nTenors, int precision = 4)
{
  // print out the results
  std::cout<< std::setw(7) << "| Delta\\T | ";
  for ( int i = 0; i < nTenors; i++ )
    std::cout << std::fixed << std::setprecision(2) << std::setw(7) << tenors[i] <<" | ";
  std::cout << std::endl;
  for ( int j = 0; j < nDeltas; j++ ) {
      std::cout<<"| " <<  std::setprecision(2) << std::setw(7) << 0.05 * (j+1) << " | ";
      for ( int i = 0; i < nTenors; i++ )
         std::cout << std::fixed << std::setprecision(precision) << std::setw(7) << v[i][j] <<" | ";
      std::cout << std::endl;
  }
}

void CalibMC(double spot, double rd, double rf, const ImpliedVol& iv, const LocalVol& lv)
{
  // generate the calibration error matrix
  std::cout << "========== MC Local Volatility Calibration ========================" << std::endl;
  const int nTenors = 15;
  const int nDeltas = 19;
  double tenors[nTenors] = { 0.02, 0.04, 0.06, 0.08, 0.16, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0};
  double bsPrices[nTenors][nDeltas];
  double strikes[nTenors][nDeltas];
  double lvPrices[nTenors][nDeltas];
  double calibErrors[nTenors][nDeltas];
  int nTPerYear = 50;// change
  int nPaths = 256*1024; // 256; //  * 1024; // *1024; // 1 million path to make sure it converges
  for ( int i = 0; i < nTenors; i++ ) {
      double T = tenors[i];
      double fwd = spot * exp((rd-rf)*T);
      vector<double> ks(nDeltas);
      for ( int j = 0; j < nDeltas; j++ ) {
         double d = 0.05 * (j+1);
         double K = StrikeFromDelta(spot, rd, rf, d, T, iv);
         strikes[i][j] = K; // just to record the strike for illustration
         // the implied volatility we can query from the implied vol surface
         double sigma = iv.Vol(T, K);
         bsPrices[i][j] = bsPricer(Call, K, T, spot, sigma, rd, rf);
         ks[j] = K;
      }
      // for mc we price a vector of trades
      auto call = [ks](double S)
                    {
                        vector<double> payoff(ks.size());
                        for ( int i = 0; i < ks.size(); i++ )
                            payoff[i] = S>ks[i]?S-ks[i]:0;
                        return payoff;
                    };
      int nT = max((int)(nTPerYear*T), nTPerYear);
      auto mcResult = mcLocalVol(call, nPaths, nT, spot, T, rd, rf, lv);
      for ( int j = 0; j < nDeltas; j++ ) {
        lvPrices[i][j] = mcResult[j].first;
        calibErrors[i][j] = (lvPrices[i][j] - bsPrices[i][j]) / spot * 10000;
	std::cout << i << ", " << j << ", error (MC, nT per year = " << nTPerYear << ", nPaths = " << nPaths << ": " << calibErrors[i][j] << ",\t stdErr(in bp): " << mcResult[j].second / spot * 10000 << std::endl; // to show progress since MC is slow
      }
  }

  // print out the results
  std::cout << "-----------------Strikes----------------------" << std::endl;
  ShowArray2(strikes,  tenors, nDeltas, nTenors);
  std::cout << "-----------------BS Price----------------------" << std::endl;
  ShowArray2(bsPrices, tenors, nDeltas, nTenors);
  std::cout << "-----------------LV Price (MC)----------------------" << std::endl;
  ShowArray2(lvPrices, tenors, nDeltas, nTenors);
  std::cout << "-----------------Calibration Error (MC, nT per year = " << nTPerYear << ", nPaths = " << nPaths << ")----------------------" << std::endl;
  ShowArray2(calibErrors, tenors, nDeltas, nTenors, 1);
}


int main(int argc, char* argv[])
{
  std::string inputFile = argc <= 1? std::string("usdcny.txt") : std::string(argv[1]);
  ifstream fin(inputFile);
  if(!fin) {
      cerr << "input file does not exist" << endl;
      return 1;
  }
  double spot, rd, rf;
  fin >> spot >> rd >> rf;
  std::cout << "-------------------Input------------------\n"
            << "spot = " << spot << std::endl
            << "rd = " << rd << std::endl
            << "rf = " << rf << std::endl
            << "t \t ATM \t MS25 \t RR25 \t MS10 \t RR10" << std::endl;
  vector< pair<double, shared_ptr<Smile> > > pillarSmiles;
  double kmin, kmax;
  while(!fin.eof()) {
      double t;
      double atmvol;
      double ms25;
      double rr25;
      double ms10;
      double rr10;
      if (fin >> t >> atmvol >> ms25 >> rr25 >> ms10 >> rr10) {
        cout << t << "\t" << atmvol << "\t" << ms25  << "\t" << rr25 << "\t" << ms10 << "\t" << rr10 << std::endl;

        vector<pair<double, double>> marks = input2Marks(spot, rd, rf, t, atmvol, ms25, rr25, ms10, rr10);
        // shared_ptr<Smile> sm_ptr(new CubicSplineSmile(marks));
        shared_ptr<Smile> sm_ptr(new ArbitrageFreeSmile(marks, spot*exp((rd-rf)*t), t));
        sm_ptr -> Construct();

        pillarSmiles.push_back( pair<double, shared_ptr<Smile>>(t, sm_ptr) );
        kmin = marks.front().first; // for plotting the charts only
        kmax = marks.back().first;
      }
  }
  ImpliedVol iv(pillarSmiles);
  LocalVol lv(iv, spot, rd, rf);
  ShowIVandLV(iv, lv, kmin, kmax); 
  CalibMC(spot, rd, rf, iv, lv);
  return 0;
}
