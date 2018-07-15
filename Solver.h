#ifndef _SOLVER_H
#define _SOLVER_H

#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>

bool RootBracketing( std::function<double(double)> f, double &a, double &b
		     , double min=std::numeric_limits<double>::min()
		     , double max=std::numeric_limits<double>::max())
{
  const int NTRY=50;
  const double FACTOR=1.6;
  if (a >= b) throw("wrong input a and b in RootBracketing");
  double f1 = f(a);
  double f2 = f(b);
  for (int j=0;j<NTRY;j++) {
    if (f1*f2 < 0.0) return true;
    if (std::abs(f1) < std::abs(f2))
      f1=f(a = std::max(a + FACTOR*(a-b), min));
    else
      f2=f(b = std::min(b + FACTOR*(b-a), max));
  }
  return false;
}

double rfbisect(std::function<double(double)> f, double a, double b, double tol)
{
  assert(a < b && f(a) * f(b) < 0);
  double c;
  while( (b-a)/2 > tol ) {
    c = (a+b) / 2.0;
    if(std::abs(f(c)) < tol)
      return c;
    else {
      if(f(a)*f(c) < 0)
	b = c;
      else
	a = c;
    }
  }
  return c;
}

double rfsecant(std::function<double(double)> f, double a, double b, double tol)
{
  double c;
  const int maxIter = 100;
  int nIter = 0;
  while( std::abs(a - b) > tol && nIter <= maxIter ) {
    c = (a * f(b) - b * f(a)) / (f(b) - f(a));
    if(std::abs(f(c)) < tol)
      return c;
    else {
      a = b;
      b = c;
    }
    nIter++;
  }
  return c;
} 

double rffalsi(std::function<double(double)> f, double a, double b, double tol)
{
  assert(a < b && f(a) * f(b) < 0);
  double c;
  while( std::abs(a - b) > tol) {
    c = (a * f(b) - b * f(a)) / (f(b) - f(a));
    if(std::abs(f(c)) < tol)
      return c;
    else {
      if(f(a)*f(c) < 0)
	b = c;
      else
	a = c;
    }
  }
  return c;
}

double rfbrent(std::function<double(double)> f, double a, double b, double tol)
{
  const int ITMAX=100;
  const double EPS=std::numeric_limits<double>::epsilon();
  
  double c=b,d,e,fa=f(a),fb=f(b); //,fc,p,q,r,s,tol1,xm;
  assert(fa * fb <= 0);
  double fc=fb, p, q, r, s, tol1, xm;
  for (int iter=0;iter<ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (std::abs(fc) < std::abs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*std::abs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (std::abs(xm) <= tol1 || fb == 0.0) return b;
    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=std::abs(p);
      double min1=3.0*xm*q-std::abs(tol1*q);
      double min2=std::abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;    
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (std::abs(d) > tol1)
      b += d;
    else                          
      b += xm>0? tol1*xm : tol1*(-xm);
    fb=f(b);
  }
  throw("Maximum number of iterations exceeded in rfbrent");
}                                                  
	
// std::function examples -- need to include "BSAnalytics.h"
// int main()
// {
  // std::function<double(double)> fun = foo;
  // std::cout << fun(5) << std::endl;

  // referring to member function
  // std::function<double(VolToBSPrice&, double)> fun2;
  // fun2 = &VolToBSPrice::operator();
  // VolToBSPrice f(100, 0.02, 90, 1.0);
  // std::cout << fun2(f, 0.04) << std::endl;

  // function partial application
  // std::function<double(double)> fun3 = std::bind(bsPricer, Call, 90, 1.0, 100, std::placeholders::_1, 0.02);
  // std::cout << fun3(0.04) << std::endl;


  // auto fun4 = [](double vol){return bsPricer(Call, 90, 1.0, 100, vol, 0.02) - 11;};
  // double a = 0.1, b = 0.5;
  // RootBracketing(fun4, a, b, 0.0001);
  // std::cout << a << " " << b << std::endl;


  // double price = bsPricer(Call, 90, 1.0, 100, 0.1, 0.02); // bs price for 10% vol
  // auto fun5 = [price](double vol){return bsPricer(Call, 90, 1.0, 100, vol, 0.02) - price;};
  // std::cout << rfbisect(fun5, a, b, 1e-6) << std::endl;

  // std::cout << rfsecant(fun5, a, b, 1e-6) << std::endl;

  // std::cout << rffalsi(fun5, a, b, 1e-6) << std::endl;

  // std::cout << rfbrent(fun5, a, b, 1e-6) << std::endl;

  // return 0;
// }

#endif
