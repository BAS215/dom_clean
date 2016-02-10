#ifndef spinorbit_
#define spinorbit_
#include "potPara.h"
#include "imaginaryForm.h"
#include <iostream>
#include <complex>

using namespace std;


/**
 *\brief spin orbit optical potential
 *
 * calculates the real and imaginary spin orbit potential
 */


class spinOrbit
{
 public:
  double Vzero;
  double V;
  double AW;
  double BW;
  double W;
  double DerV;
  double R;
  double a;
  double Efermi;
  double Ecm;
  potPara Real;
  potPara Disp;
  potPara DerDisp;
  potPara Imag;
  imaginaryForm Form;
  void load(double,double,double,double,double,double);
  void load(double V0, double W0, double R0, double a0);
  void setEnergy(double);
  double RealPotential(double);
  double ImaginaryPotential(double);
  double DispersiveCorrection(double r);
  double DerDispersiveCorrection(double r);
  complex<double>potential(double r);
  complex<double>potentialE(double r, double Ecm);
  double setAM(int l, double j);
  double LdotSigma;
};

#endif
