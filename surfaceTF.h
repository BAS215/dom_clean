#ifndef surfaceTF_
#define surfaceTF_

#include "potPara.h"
#include "twoFermi.h"
#include <iostream>
#include <complex>

using namespace std;

/**
 *\brief imaginary surface potential and dispersive correction
 *
 *this class deals with all aspects of the surface imaginary potential, its
 *dispersive correction and contributions to effective mass and occupation
 *probabilities. the magnitude of the imaginary potential is parametrized as
 *\f$ W(E) = A \frac{\|X_{1}\|^{m}}{\|X_{1}\|^{m}+B^m} \frac{1}{1+\exp\left( \frac{X_{1}-C}{D}\right)} \frac{1}{1+\exp-\left(\frac{X_{1}+C}{D}\right)} \f$ 
 * \f$ X_{1} = E-E_{Fermi} \f$
 */


class surfaceTF
{
 private:
  static double const konst;
 public:

  surfaceTF();
  surfaceTF(int typeN0);
  complex<double> potentialE(double r1, double r2, double dr, double Ecm);
  complex<double> potentialE(double r, double Ecm);

  complex<double> potential(double r1, double r2, double dr);
  complex<double> potential(double r);

  double H(double dr);
  complex<double> U(double r1, double r2);
  double dU(double r1, double r2);

  void load(double,double,double,double,double,double,double,double,double,int);
  void load(double strength, double R0, double a0, double beta0);
  void setEnergy(double);
  double imaginaryE( double E ); 
  double dispersiveE( double E ); 
  double derDispersiveE( double E ); 
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerDispersiveCorrection(double r);


  double CentralImaginaryPotential(double);
  double CentralDeltaRealPotential(double);
  double getMaxW();
  void setTypeN(int typeN0);

  potPara form;  //wood saxon form



  twoFermi TF; // energy dependence of imaginary potential

  double A;
  double B;
  double C;
  double D;
  double Wstart;
  double Efermi; //!< Fermi energy im MeV
  double Ecm;

  double R; //!> radius of potential
  double a; //!< diffuseness of potential
  double beta; //!< nonlocality parameter

  double V; //!< real potential for a given energy
  double W; //!< imaginary potential for a given energy
  double derV; //!< magnitude of E derivative of V

  int typeN; //!< which nonlocatility U function to use
  int EnergyType; //!< Which Dispersive correction to use

};

#endif 
