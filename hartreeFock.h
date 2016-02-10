#ifndef hartreefock_
#define hartreefock_

#include <stdio.h>
#include <stdlib.h>
#include "potPara.h"
#include <cmath>
#include <iostream>
#include "Gauss_q.h"

using namespace std;

/**
 *\brief give the energy-dependent Hartree-Fock potential
 *
 *Deals with all aspects of the Hartree-Fock potential - This is the
 * main real potential without the dispersive corrections. The energy
 * dependece is to take into account nonlocality
 */
class hartreeFock
{
 private:
  static double const konst; 

public:

  hartreeFock();
  hartreeFock(int typeN0);

  void load(double Vvol, double Vsur, double R0, double a0, double, double, double beta);

  void load(double Vvol, double Vsur, double R0, double a0,  double , double , double beta, double A, double beta20);

  void load(double Vvol, double Vsur, double R0, double a0, double beta, double V_wine1, double R_wine1, double rho_wine1);

  void load(double Vvol, double Vsur, double R0v, double a0v, double R0s, double a0s, double beta, double A, double beta20 , double V_wine1, double R_wine1, double rho_wine1);

  double potential(double r);
  double potential0(double r);
  double potential1(double r);
  double potentialSurface(double r);
  double potentialEquivLocal(double r);
  double potential (double r1, double r2, double dr);
  double potential (double*R1, double *R2);
  void setEquivLocalDepth(double VvolE);
  void setTypeN( int typeN0);
  double U(double r1, double r2);
  double U0(double r1, double r2);
  double U1(double r1,double r2);

  double BWU_local(double);
  double BWU(double r1, double r2);

  double U_Angular_Integration(double, double ,double, int);

  double H(double dr, double beta);

  potPara volume;  //!< volume type form factor
  potPara surface; //!< surface type form factor
  potPara bottleW;

  double Rv;  //!< radius of potential in fm
  double av;  //!< diffuseness of potential in fm
  double Rs;  //!< radius of potential in fm
  double as;  //!< diffuseness of potential in fm
  double Vvol; //!< depth of volume component in MeV
  double beta0; //!< nonlocality component
  double Vsur; //!< strength of surface component in MeV

  double A;  //!< magnitude of second componet of nonlocality 
  double beta1; //!< nonlocality length of second component

  double VvolE; //!< depth of equivalent local potential
   
  double V_wine;
  double R_wine;
  double rho_wine;
  int typeN; //!< which nonlocatility U function to use

};

#endif
