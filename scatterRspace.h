#include "pot.h"
#include "numerov.h"
#include "waves.h"
#include "legendre.h"
#include "compound.h"
#include <string>

using namespace std;

/**
 *!\brief   nonlocal calculation of scattering in r space
 *
 * Solves the (non local) Schrodinger Equation for position energy.
 * from this gets the S matrix and then one can determined the elastic
 * scattering angular distribution, the reaction cross section, etc.
 * used the iterative technique of N. Michel (EPJ ) to solve the Schrodinger 
 * equation in position space 
 */



class scatterRspace : public numerov, public compound
{
 private:
  static double const pi;

  static double const kconstant;
  static double const e2;

  double coulDis;
  void nonLocalArray1(int l);
  void nonLocalArray0(int l);

 public:
  static double const m0;
  double konst;
  double Z; //!< atomic number of target
  double Zp; //!<atomic number of nucleon
  double A; //!< mass number of target
  double Kwave; //!< asymptotic wavelength
  double Kwave2; //!< square of Kwave
  double gammaRel; //!< relativistic gamma factor
  double gamma; //!<Sommerfeld parameter
  double mu;//!< reduced mass
  double muhbar; //<! 2*mu/hbar**2
  double Ecm; //!< center-of-mass energy in MeV
  double Elab; //!< lab energy in MeV

  void init(double Z0, double Zp0, double A0);
  scatterRspace(double Z0, double Zp0, double A0, pot*Pot, string * title);
  ~scatterRspace();
  void setPot(pot*Pot0);
  pot * Pot;
  int getSmatrix(double Ecm, double Elab);
  double DifferentialXsection(double);
  double AnalyzePower; // Analysing power, calculated by DifferentialXsection
  double SpinRotation; // spin rotation parameter Q, calculated by Differential
  double AbsorptionXsection(); //absorption cross section
  double TotXsection(); //total cross section (for neutrons)
  double ElasticXsection(); // total elastic cross section (for neutrons)
  double Rutherford(double);
  
  void integrateWaveFunction(int istart, int l);

  int lMax;  //!< maximum orbital am that can be considerd.
  int llMax; //!< array dimension for  orbital stuff
  int lStop; //!< maxiumin orbital AM for which the S matrix is calculated
  double **phaseShift; // array of dimension [lMax+1][2] with phaseshifts;
  double **S2; //array of dimension [lMax+1][2];
  complex<double> **S; //array of dimension [lMax+1][2];
  double *Sigma;//array of dimension [lMax+1] containing Coulomb phase shifts

  int Nnumerov; //!< number of point in r array for Numerov integration
  double rStart; //!<radius to start integration (fm)
  double rStop; // !<radius to stop interation (fm) 


  double * SigmaAb; //!<absorption cross section for each j
  double * SigmaAb_pos; // !<for postive parity
  double * SigmaAb_neg; // !<for negative parity

  double TransCoef(int, double);
  void newEnergy(double Ecm, double Elab);
  double mstar;//!< effective mass relative to nucleon mass
};

