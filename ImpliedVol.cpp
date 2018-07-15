#include "ImpliedVol.h"
#include "LocalVol.h"
#include "BSAnalytics.h"
//#include "Solver.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <random>
#include <cassert>
#include <cmath>

#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

using namespace boost::numeric;
using namespace std;

// choose exact integral type
//#ifdef CGAL_USE_GMP
//#include <CGAL/Gmpz.h>
//typedef CGAL::Gmpz ET;
//#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
//#endif

void CubicSplineSmile::Construct()
{
  // if we inject two points to smooth the wings
  pair<double, double> leftNode (marks[0].first / 2.0, 
				 marks[0].second - marks[0].first / 4.0 * (marks[1].second-marks[0].second) / (marks[1].first-marks[0].first) );
  pair<double, double> rightNode (2*marks[marks.size()-1].first - marks[marks.size()-2].first, 
				  marks[marks.size()-1].second + (marks[marks.size()-1].second-marks[marks.size()-2].second) / 2.0 );
  marks.insert(marks.begin(), leftNode);
  marks.push_back(rightNode);
  
  int N = marks.size();
  // end y' are zero, flat extrapolation
  double yp1 = 0.0; // 01;
  double ypn = 0.0; // 01;
  y2.resize(N);
  vector<double> u(N-1);

  y2[0] = -0.5;
  u[0]=(3.0/(marks[1].first-marks[0].first)) *
    ((marks[1].second-marks[0].second) / (marks[1].first-marks[0].first) - yp1);

  for(int i = 1; i < N-1; i++) {
    double sig=(marks[i].first-marks[i-1].first)/(marks[i+1].first-marks[i-1].first);
    double p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(marks[i+1].second-marks[i].second)/(marks[i+1].first-marks[i].first)
      - (marks[i].second-marks[i-1].second)/(marks[i].first-marks[i-1].first);
    u[i]=(6.0*u[i]/(marks[i+1].first-marks[i-1].first)-sig*u[i-1])/p;
  }

  double qn=0.5;
  double un=(3.0/(marks[N-1].first-marks[N-2].first)) *
    (ypn-(marks[N-1].second-marks[N-2].second)/(marks[N-1].first-marks[N-2].first));

  y2[N-1]=(un-qn*u[N-2])/(qn*y2[N-2]+1.0);
  for (int i=N-2;i>=0;i--) {
    y2[i]=y2[i]*y2[i+1]+u[i];
  }
}

double CubicSplineSmile::Vol(double strike) const
{
  int i;
  // we use trivial search, but can consider binary search for better performance
  for (i = 0; i < marks.size(); i++ ) 
    if (strike < marks[i].first )
      break; // i stores the index of the right end of the bracket
  // std::cout << " i = " << i << std::endl;
  // for (i = 0; i < marks.size(); i++ ) 
  // std::cout << marks[i].first << ", " << marks[i].second << std::endl;

  // extrapolation
  if (i == 0)
    return marks[i].second; //  - (strike - marks[i].first) * 0.001 ;
  if (i == marks.size() )
    return marks[i-1].second; //  + (strike - marks[i-1].first) * 0.001 ;

  // interpolate
  double h = marks[i].first - marks[i-1].first;
  double a = (marks[i].first - strike) / h;
  double b = 1 - a;
  double c = (a*a*a - a) * h * h / 6.0;
  double d = (b*b*b - b) * h * h / 6.0;
  return a*marks[i-1].second + b*marks[i].second + c*y2[i-1] + d*y2[i];
}

void ArbitrageFreeSmile::Construct()
{
  // program and solution types
typedef CGAL::Nonnegative_quadratic_program_from_iterators
<double**,                                                // for A
 double*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>,    // for r
 double**,                                                // for D
 double*>                                                 // for c 
 Program;
 typedef CGAL::Quadratic_program_solution<ET> Solution;

  const int N = 50;   // number of sample points
  const int M = marks.size();
  // calculate the length of segment
  double sigmaATM = marks[(M + 1) / 2].second; // sigmaATM is the volatility of the 
                                               // middle point of the input marks
  double sigmaSqrtT = sigmaATM * std::sqrt(T);
  double tmp1 = - sigmaSqrtT * sigmaSqrtT / 2;
  double tmp2 = 5 * sigmaSqrtT;
  const double u = F * (exp(tmp1 + tmp2) - exp(tmp1 - tmp2)) / (N - 1);
  double CallPrice[5]={0};   // opition price of 5 marks
  // calculate the strike price and length of segment
  K.resize(N);
  for (int i = 0; i < N; i++) 
    {
		K[i] = F * exp(tmp1 - tmp2) + i * u;
    }
  
  for(int i = 0; i < marks.size(); i++)
	{
		CallPrice[i]= bsPricer(Call, marks[i].first, T, F, marks[i].second);	
	}
	
	double Ac[N-2][N+4], Ap[N-2][N+4];  // eliminate column 1 and N and move c1 and cN to RHS to achieve constraint 5
	for(int i = 0; i < N-2; i++)
	{
		for(int j = 0; j < N+4; j++)
		{
			Ac[i][j] = 0;
			Ap[i][j] = 0;
		}
	}
	
	double b[N+4] = {0};    
	b[0] = - F + K[0];
	b[48] = 1.0 / u;  // constraint 4 
	Ac[0][0] = -2; 
	Ac[0][1] = 1;
	Ac[47][46] = 1;
	Ac[47][47] = -2;
	Ap[0][0] = -2.0 / 3.0 * u * u; 
	Ap[0][1] = -1.0 / 6.0 * u * u;
	Ap[47][46] = -1.0 / 6.0 * u * u;
	Ap[47][47] = -2.0 / 3.0 * u * u; 
	Ap[0][48] = 1;
	Ap[47][48] = 1;
	for(int i = 1; i <= 46; i++)
	{
		Ac[i][i-1] = 1;
		Ac[i][i] = -2;
		Ac[i][i+1] = 1;
		Ap[i][i-1] = -1.0 / 6.0 * u * u;
		Ap[i][i] = -2.0 / 3.0 * u * u;
		Ap[i][i+1] = -1.0 / 6.0 * u * u;
		Ap[i][48] = 1;
	}

// constraint 2
double a_[5] = {0};
double b_[5] = {0};
for(int i = 0; i < marks.size(); i++)
{
	int j;
	for(j = 1; j < N; j++)
	{
		if(marks[i].first < K[j])
			break;
	}
	if(j == 1) 
	{
		a_[i] = (K[j] - marks[i].first) / u;
		b_[i]=1 - a_[i];
		Ac[j-1][48+i+1] = b_[i];
		Ap[j-1][48+i+1] = (b_[i] * b_[i] * b_[i] - b_[i]) / 6.0 * u * u;
		b[48+i+1] = CallPrice[i] - a_[i] * (F - K[0]);
	}
	if(j == N-1)
	{
		a_[i] = (K[j] - marks[i].first) / u;
		b_[i] = 1 - a_[i];
		Ac[j-2][48+i+1] = a_[i];
		Ap[j-2][48+i+1] = (a_[i] * a_[i] * a_[i] - a_[i]) / 6.0 * u * u;
		b[48+i+1] = CallPrice[i];
	}
	else
	{
		a_[i] = (K[j] - marks[i].first) / u;
		b_[i] = 1 - a_[i];
		Ac[j-2][48+i+1] = a_[i];
		Ac[j-1][48+i+1] = b_[i];
		Ap[j-2][48+i+1] = (a_[i] * a_[i] * a_[i] - a_[i]) / 6.0 * u * u;
		Ap[j-1][48+i+1] = (b_[i] * b_[i] * b_[i] - b_[i]) / 6.0 * u * u;
		b[48+i+1] = CallPrice[i];
	}
}

double *A[]={Ac[0],Ac[1],Ac[2],Ac[3],Ac[4],Ac[5],Ac[6],Ac[7],Ac[8],Ac[9],Ac[10],Ac[11],Ac[12],Ac[13],Ac[14],Ac[15],Ac[16],Ac[17],Ac[18],Ac[19],Ac[20],Ac[21],Ac[22],Ac[23],Ac[24],
          Ac[25],Ac[26],Ac[27],Ac[28],Ac[29],Ac[30],Ac[31],Ac[32],Ac[33],Ac[34],Ac[35],Ac[36],Ac[37],Ac[38],Ac[39],Ac[40],Ac[41],Ac[42],Ac[43],Ac[44],Ac[45],Ac[46],Ac[47],
          Ap[0],Ap[1],Ap[2],Ap[3],Ap[4],Ap[5],Ap[6],Ap[7],Ap[8],Ap[9],Ap[10],Ap[11],Ap[12],Ap[13],Ap[14],Ap[15],Ap[16],Ap[17],Ap[18],Ap[19],Ap[20],Ap[21],Ap[22],Ap[23],Ap[24],
		  Ap[25],Ap[26],Ap[27],Ap[28],Ap[29],Ap[30],Ap[31],Ap[32],Ap[33],Ap[34],Ap[35],Ap[36],Ap[37],Ap[38],Ap[39],Ap[40],Ap[41],Ap[42],Ap[43],Ap[44],Ap[45],Ap[46],Ap[47]};

CGAL::Const_oneset_iterator<CGAL::Comparison_result> 
        r(    CGAL::EQUAL);                 // constraints are "="	
double h[2*N-4][2*N-4];   // 96*96
for(int i=0; i < 2*N-4; i++)
{
	for(int j=0; j < 2*N-4; j++)
		h[i][j]=0;	
}

h[48][48] = 4.0 / 3.0 * u * u; 
h[48][49] = 1.0 / 3.0 * u * u;
h[95][94] = 1.0 / 3.0 * u * u; 
h[95][95] = 4.0 / 3.0 * u * u;

for(int i=49;i<=94;++i)
{
		h[i][i-1] = 1.0 / 3.0 * u * u;
		h[i][i] = 4.0 / 3.0 * u * u;
		h[i][i+1] = 1.0 / 3.0 * u * u;
}
double* D[] = {h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],h[11],h[12],h[13],h[14],h[15],h[16],h[17],h[18],h[19],h[20],h[21],h[22],h[23],h[24],h[25],
          h[26],h[27],h[28],h[29],h[30],h[31],h[32],h[33],h[34],h[35],h[36],h[37],h[38],h[39],h[40],h[41],h[42],h[43],h[44],h[45],h[46],h[47],h[48],h[49],h[50],h[51],
		  h[52],h[53],h[54],h[55],h[56],h[57],h[58],h[59],h[60],h[61],h[62],h[63],h[64],h[65],h[66],h[67],h[68],h[69],h[70],h[71],h[72],h[73],h[74],h[75],h[76],h[77],
		  h[78],h[79],h[80],h[81],h[82],h[83],h[84],h[85],h[86],h[87],h[88],h[89],h[90],h[91],h[92],h[93],h[94],h[95]};
double c[96] = {0};
double c0 = 0;
  
  // now construct the quadratic program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program qp (96, 54, A, b, r, D, c, c0);
	
  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());
  
  // output solution
   std::vector<double> sol;
  for (  Solution::Variable_value_iterator iter = s.variable_values_begin(); iter != s.variable_values_end(); iter++) 
  {
    sol.push_back(CGAL::to_double(*iter));
  }

  // output solution to class member
  C.resize(K.size(), .0);
  y2.resize(K.size(), .0);
  C[0]=F-K[0];C[49]=0;y2[0]=0;y2[49]=0;
  for(int i=0;i<=47;++i)
  {
	C[i+1]=sol[i];
	y2[i+1]=sol[48+i];
  }

  // output the results to std::cout
  std::cout << std::endl;
  std::cout << "Status: ";
  switch(s.status()) {
    case 1:
      std::cout << "INFEASIBLE" << std::endl;
      break;
    case 2:
      std::cout << "UNBOUNDED" << std::endl;
      break;
    case 3:
      std::cout << "OPTIMAL" << std::endl;
      break;
    default:
      break;
  }
  std::cout << std::endl;

  if (s.status() == 3) {
    std::cout << "Variables:" << std::endl;
    std::cout << "C: ";
    for (int i = 0; i < (int)C.size(); ++i) {
      std::cout << C[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "y2: ";
    for (int i = 0; i < (int)y2.size(); ++i) {
      std::cout << y2[i] << " ";
    }
    std::cout << std::endl;
  }
}

double ArbitrageFreeSmile::Vol(double strike) const
{
  assert(K.size() == C.size());

  double price(0.0);
  int l = SearchStrike(strike);
  if (l == 0) {
    price = C[0];
  } else if (l == K.size()) {
    price = C[C.size() - 1];
  } else {
    double h = C[l] - C[l - 1];
    double a = (K[l] - strike) / h;
    double b = 1 - a;
    double c = (a * a * a - a) * h * h / 6.0;
    double d = (b * b * b - b) * h * h / 6.0;
    price = a * C[l - 1] + b * C[l] + c * y2[l - 1] + d * y2[l];
  }

  if (price > bsPricer(Call, strike, T, F, 0)) {
    double sigma = bsImpliedVol(Call, strike, T, F, price);
    return sigma;
  } 
  else 
  {
    return 0.0;
  }
}

int ArbitrageFreeSmile::SearchStrike(double strike) const
{
  if (K.empty()) return -1;
  // we use trivial search, but can consider binary search for better performance
  for (int i = 0; i < (int)K.size(); ++i) {
    if (strike < K[i]) {
      return i;
    }
  }
  // strike is larger than all of the K
  return K.size();
}

ImpliedVol::ImpliedVol(const vector<pair<double, shared_ptr<Smile>>>& _pillarSmiles)
  : pillarSmiles(_pillarSmiles)
{
}

double ImpliedVol::Vol(double t, double k) const
{
  // we use trivial search, but can consider binary search for better performance
  int i;
  for (i = 0; i < pillarSmiles.size(); i++ ) 
    if (t < pillarSmiles[i].first )
      break; // i stores the index of the right end of the bracket

  if (i == 0)
    return pillarSmiles[0].second->Vol(k);
  if (i == pillarSmiles.size())
    return pillarSmiles[i-1].second->Vol(k);
  double t1 = pillarSmiles[i-1].first;
  double t2 = pillarSmiles[i].first;
  double a = (t2 - t) / (t2 - t1);
  double b = 1 - a;
  double v1 = pillarSmiles[i-1].second->Vol(k);
  double v2 = pillarSmiles[i].second->Vol(k);  
  return std::sqrt( (a*v1*v1*t1 + b*v2*v2*t2)/t );
}
