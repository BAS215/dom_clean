#ifndef _numerov
#define _numerov
#include <complex>
#include <cmath>
#include <iostream>


using namespace std;
class numerov
{
 private:
  double xstart;
  double xend;


  double fact;
  bool start;

 public:
  int N;
  double dx;
  double *x;
  complex<double> *y;
  complex<double> *f;
  complex<double> *source;


  numerov();
  void initNumerov(double xstart, double xend, int N);
  virtual ~numerov();
  void solveNumerov0(complex<double> y0, complex<double> dydx0);
  void solveNumerov1(complex<double> y0, complex<double> y1);
  void solveNumerov(int istart);
};
#endif
