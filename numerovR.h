#ifndef _numerovR
#define _numerovR
#include <cmath>
#include <iostream>


using namespace std;
class numerovR
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
  double *y;
  double *f;
  double *source;
  double *vloc;

  numerovR();
  void initNumerovR(double xstart, double xend, int N);
  virtual ~numerovR();
  void solveNumerovR(double y0, double dydx0);

  bool generateVloc;
};
#endif
