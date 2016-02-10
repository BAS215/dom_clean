#include "minimize1D.h"
#include <cmath>
#include <iostream>
using namespace std;

//parabolic interpolation to find minimum
//translated from "Numerical recipes in Fortran 77",
//Press, Teukolsky, Vetterling, Flannery, 2nd edition Page 397
double minimize1D::Brent(double*y,const double& tolerance)
{
  const double Zeps = 1.e-10;
  const double gold = .3819660;
  const int itMax = 100;
  int Golden; // flag for golden rule
  //a and b must be in ascending order
  double a,b;
  if (y[0] > y[2]) 
    {
      a = y[2];
      b = y[0];
    }
  else 
    {
      b = y[2];
      a = y[0];
    }

  //initialize
  double v = y[1];
  double w = v;
  double x = v;
  double e = 0.;  // this will be the distance moved on the step before last
  double fx = funct1D(x);
  double fv = fx;
  double fw = fx;
  double d=0.;
  double u,fu;
  double fuu = -1.;
  for (int iter=0;iter<itMax;iter++)
    {
      double xm = (a+b)/2.;
      double tol1 = tolerance*abs(x)+Zeps;
      double tol2 = 2.*tol1;
      //test for completion
      if (abs(x-xm) <= (tol2-0.5*(b-a))) 
	{
	  y[0] = x;
          return fx;
	}
      Golden = 1;
      if (abs(e) > tol1)
	{
	  double r = (x-w)*(fx-fv);
	  double q = (x-v)*(fx-fw);
	  double p = (x-v)*q-(x-w)*r;
          q = 2.*(q-r);
	  if (q > 0.)p = -p;
	  q = abs(q);
	  double etemp = e;
	  e = d;
          // checl for acceptability of parabolic fit
	  if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) 
	    Golden = 1; // fit ok
	  else // fit not ok
	    {
	      Golden = 0;  //skip over the golden section step
	      d = p/q;
	      u = x + d;
              if (u-a < tol2 || b-u < tol2) 
		{
		  if (xm-x > 0.)d = tol1;
		  else d = -tol1;
		}
	    }

	}
      //we arrive here at the Golden section step, which we take into the 
      //larger of the two segments
      if (Golden)
	{
	  if (x > xm)e = a - x;
	  else e = b - x;
	  d = gold*e; // take the golden section
	}
      //at this point d is calculated either from the parabolic fit or the 
      //golden section
      if (abs(d) >= tol1) u = x + d;
      else 
	{
	  if (d > 0.) u = x + tol1;
	  else u = x - tol1;
	}
      fu = funct1D(u); // function evalution per step
      //housekeeping

      if (fu == fx && fuu == fx) 
	{
	  cout << "breaking out" << endl;
         break;
	}

      if (fu <= fx) 
	{
          if (u >= x) a = x;
          else b = x;
          v = w;
          fv = fw;
          w = x;
          fw = fx;
          x = u;
          fuu = fx;
          fx = fu;

	}
      else
	{
	  if ( u < x) a = u;
	  else b = u;
	  if (fu <= fw || w == x)
	    {
	      v = w;
	      fv = fw;
	      w = u;
	      fw = fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    {
	      v = u;
	      fv = fu;
	    }

	}
    }
  cout << "exceed maximum iteractions" << endl;
  y[0] = x;
  return fx;
}
//************************************************************************* 
// brackets the minimum
// given initial points x[0] and x[2], this returns
// three (possibly different) points x[0],x[1],x[2]
// such that f[1]=funct(x[1]) is smaller than either
// f[0]=funct(x[0]) or f[2]=funct(x[2]). Thus a minimum must 
// be between x[0] and x[2]. 
// translated from "Numerical recipes in Fortran 77" Press, Teukolsky, 
//Vetterling, and Flannery, 2nd edition Page 393
void minimize1D::bracket(double*x,double*f)
{
  const double Rlimit = 100.;
  const double tiny = 1.e-20;
  const double ratio = 1.618034;
  f[0] = funct1D(x[0]);
  f[1] = funct1D(x[1]);
  double fu;

  // want f[1] < f[0], if not, swap them
  if( f[1] > f[0])
    {
      double temp = x[0];
      x[0] = x[1];
      x[1] = temp;
      temp = f[0];
      f[0] = f[1];
      f[1] = temp; 
    }

  //first guess for x[2]
  x[2] = x[1] + ratio*(x[1]-x[0]);
  f[2] = funct1D(x[2]);
  while (f[1] >= f[2])
    {
      //use parabolic extrapolation to find point u
      double r = (x[1]-x[0])*(f[1]-f[2]);
      double q = (x[1]-x[2])*(f[1]-f[0]);
      double denon = q-r;
      // make sure we don't divide by zero
      if (abs(denon) < tiny)
	{
	  if (denon < 0) denon = -tiny;
	  else denon = tiny;
	}
      double u = x[1] - (( x[1]-x[2])*q-(x[1]-x[0])*r)/(2.*denon);
      // set a limit that we shouldn't go futher than
      double uLimit = x[1] + Rlimit*(x[2]-x[1]);
      if ((x[1]-u)*(u-x[2]) > 0.) // u is  between x[1] and x[2], try it
	{
	  fu = funct1D(u);
          if (fu < f[2])  //we have found a minimum between x[1] and x[2]
	    {
	      x[0] = x[1];
              f[0] = f[1];
              x[1] = u;
              f[1] = fu;
              //x[2] unchanged
	      return;
	    }
	  else if (fu > f[1]) //got a minimum between x[0] and u
	    {
	      x[2] = u;
	      f[2] = fu;
              return;
	    }
          // parabolic fit didn't work, use default magnification
          u = x[2] + ratio*(x[2]-x[1]); 
          fu = funct1D(u);
	}
      else if ((x[2]-u)*(u-uLimit) > 0.) //parabolic fit is between x[2]
	                                 //and the procribed limit 
	{
	  fu = funct1D(u);
          if (fu < f[2])
	    {
              x[1] = x[2];
	      x[2] = u;
	      u = x[2] + ratio*(x[2]-x[1]);
	      f[1] = f[2];
	      f[2] = fu;
	      fu = funct1D(u);
	    }
	} 
      else if ((u-uLimit)*(uLimit-x[2]) >= 0.) // limit u to the 
                                               //maximum allowed value
	{
	  u = uLimit;
	  fu = funct1D(u);
	}
      else     //reject parabolic u, use default magnification
	{
	  u = x[2] + ratio*(x[2]-x[1]);
	  fu = funct1D(u);
	}
      //eliminate oldest point and continue
      x[0] = x[1];
      x[1] = x[2];
      x[2] = u;
      f[0] = f[1];
      f[1] = f[2];
      f[2] = fu;
    }
}

