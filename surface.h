#ifndef surface_
#define surface_

#include "potPara.h"
#include "imaginaryForm.h"
#include <iostream>

using namespace std;

/**
 *\brief imaginary surface potential and dispersive correction
 *
 *this class deals with all aspects of the surface imaginary potential, its
 *dispersive correction and contributions to effective mass and occupation
 *probabilities. the magnitude of the imaginary potenteial is parametrized as
 *\f$ W(E) = A\Theta(X_{1})\frac{X_{1}^{m1}}{X_{1}^{m1}+B_{1}^m1} \exp \left(C X_{1}\right) - A\Theta(X_{2}) \frac{X_{2}^{m2}}{X_{2}^{m2}+B_{2}^m2} 
\exp\left(C X_{2} - C (E_{p}^2-E_{p}^1)\right)\f$ where 
 * \f$ X_{1} = |E-E_{Fermi}| - E_{p}^{1} \f$ and 
 * \f$ X_{2} = |E-E_{Fermi}| - E_{p}^{2} \f$ and 
 * \f$ E_{p}^2 = B_{1} + \Delta E_{p} \f$
 */


class surface
{
 public:
  void load(double,double,double,double,double,double,double,double,double,
            int,int);
  void load(double strength, double R0, double a0);
  void SetEnergy(double);
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerivativeDispersive(double);
  double DerivativeParticleHole(double);

  double CentralImaginaryPotential(double);
  double CentralDeltaRealPotential(double);
  double CentralDerDeltaRealPotential(double);

  potPara Imaginary;  //imaginary potential 
  potPara Dispersive; //dispersive coorection to real potential
  potPara DerDispersive; //derivative of dispersive need for effective mass
  potPara ParticleHole; // derivative of particle of hole contribution to 
                   // dispersive correction needed for occupation probabilities
  imaginaryForm form1; // energy dependence of imaginary potential
  imaginaryForm form2;

  double A1;
  double B1;
  double Ep1;
  double C;
  double DeltaEp;
  double A2;
  double B2;
  double Ep2;
  double Efermi;
  double Ecm;

  double R;
  double a;
  int m1;
  int m2;

};

#endif 
