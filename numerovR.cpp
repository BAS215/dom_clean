#include "numerovR.h"


numerovR::numerovR()
{
  start = false;
}

void numerovR::initNumerovR(double xstart0, double xend0, int N0)
{
  if (start) return;
  xstart = xstart0;
  xend = xend0;
  N = N0;
  dx = (xend-xstart)/(double)N;

  fact = pow(dx,2)/12.;
  x = new double [N+2];
  y = new double [N+2];
  f = new double [N+2];
  vloc = new double [N+2];
  source = new double [N+2];


  for (int i=0;i<N+2;i++)
    {
      x[i] = xstart + (double)i*dx;
    }
  start = true;
}

//*************************************************************
numerovR::~numerovR()
{
  if (start)
    {
     delete [] x;
     delete [] y;
     delete [] f;
     delete [] source;
     delete [] vloc;

    }
}
//**************************************************************
void numerovR::solveNumerovR(double y0, double dydx0)
{



  y[0] = y0;

  y[1] = y[0]*(1.-f[2]*fact/2.) + dx*dydx0*(1.-f[2]*fact);
  y[1] += fact/2.*(7.*f[0]*y[0]+ 6.*source[1]-source[2])
    - pow(dx,4)/36.*f[2]*(f[0]*y[0]+2.*source[1]);
  y[1] /= 1. - f[1]*pow(dx,2)/4. + f[1]*f[2]*pow(dx,4)/18.;



  for (int i=2;i<N+2;i++)
    {
      y[i] = 2.*y[i-1] - y[i-2];
      y[i] += fact*(source[i]+10.*f[i-1]*y[i-1] + f[i-2]*y[i-2]);
      y[i] /= (1.-fact*f[i]);
    }

  // for (int i=0;i<N+2;i++) cout << x[i] << " " << y[i].real() << " " << y[i].imag() << endl;

}


