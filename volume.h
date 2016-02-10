


#ifndef volume_
#define volume_
#include "potPara.h"
#include "imaginaryForm.h"
#include "surVolume.h"
#include <complex>

/**
 *\brief imaginay volume potential and dispersive potential
 *
 *This class deals with all aspects of the volume imaginary potential, its
 *dispersive correction and contributions to effetive mass and occupation
 *probabilities
 *
 *there is a main volume component (with constant geometry, i.e. r 
 *dependence), plus a surface correction with acts to approximate a change of
 *radius with energy. At Ecm = Ef the radius is Rzero + deltaR,
 * at Ecm = (+-) infinity the radius is Rzero. 
 */


class volume
{
 private:
  static double const konst;
 public:
  complex<double> potential(double r);
  complex<double> potential(double r1, double r2, double dr);
  complex<double> potentialE(double r, double Ecm);
  complex<double> potentialE(double r1, double r2, double dr, double Ecm);
  complex<double> U(double r1, double r2);
  complex<double> Un();
  double dUn();
  double dU(double r1, double r2); 
  double H(double dr);
  void setTypeN(int typeN0);
  volume();
  volume(int typeN0);


  void load(double,double,double,double,double,double,double,double,
	    int,int,double,double,double,int);

  void load(double,double,double,double,double,double,double,double,
	    int,int,double,double,double,int,double,double,double);
  
  void load(double strength, double R0, double a0,double beta0);
  void setEnergy(double);
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerDispersiveCorrection(double);
  double DerivativeParticleHole(double);

  //hartreeFock hart;
  potPara form; //!< woodSaxon form factor
 
  imaginaryForm volumeF; //!< main volume contribution with fixed radius
  //imaginaryForm surface; // !<surface correction to give change of radius
  surVolume surface; //!, surface correction to give change of radius
   double deltaR; // !,energy dependence of radius
   double Rzero; // !<radius at E= Efermi
   double expR; //!<exponential decrease of radius parameter
   double a;  //!<diffuseness
   double A;  
   double B;
   double Ep; //!< energy gap either side of  fermi energy where W==0 
   double Efermi;  //!<Fermi energy
   int Asy;  //!<switch for energy-asymmetry contribution
   int m;
   double alpha; //!<paramter for energy asymmetry term
   double Ea;  //!<parameter for energy asymmetry term
   double Ecm; //!<center-of-mass Energy in MeV
   double beta; //!< nonlocality parameter 
   int typeN; //!< method to calculate U function
   int Energytype;
   double Vvol;  //!< magnitude of the real dipservive correctio
   double Wvol;  //!< magnitude of the imaginary potential
   double Vsur;  //!< magnitude of the surface correction to the real dispersive
   double Wsur; //!< magnitude of the surface correction to the imaginary potential
   double derVvol; //!< magnitude of derivative of dispersive correction
   double derVsur; //!< magnitude of surface correction to derivative of dispersive correction
   double WvolE; //!< equivalent effective nonlocal potential
   void setEquivLocalDepth(double WvolE0);
   double potentialEquivLocal(double r);
};
#endif
