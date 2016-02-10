#include "minimize1D.h"


/**
 *\brief base class for function minimization in multi-dimension
 *
 *base class for function minimization in multi-dimension
 *it itself uses the 1D minimization base class
 */

class minimizeND: public minimize1D
{
public:
  minimizeND(){}
  minimizeND(int);
  ~minimizeND();
  int ND; // dimension of function
  virtual double functND(double*)=0; // multidimensional function must be 
                                     //provided
  double* pcom;
  double* xicom;
  double funct1D(double);
  double linmin(double*,double*);
  double powell(double*,double**,double const &); // minimization function
};
