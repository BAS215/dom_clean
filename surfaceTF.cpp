#include "surfaceTF.h"

double const surfaceTF::konst  = 1./pow(acos(-1.),3./2.);

/**
 * loads in all the paramters defining the radial and energ dependence of the
 * surface imaginary potential
\param R0 is the radius parameter in fm
\param a0 is the diffuseness parameter in fm
\param A0 defines the magnitude of the potential
\param B0 is an energy-dependent parameter 
\param C0 is an energy-dependent parameter
\param D0 is an energy dependent parameter
\param Wstart0 - energy from the Fermi energy where potential starts (MeV)
\param Efermi0 - Fermi energy (MeV)
\param beta0 - nonlocality length in fm 
\param EnergyType - 0 for E < Efermi, 1 for E > Efermi, 2 for both regions
 */

void surfaceTF::load(double R0, double a0, double A0, double B0, double C0,
                     double D0, double Wstart0, double Efermi0,
                     double beta0, int EnergyType0 ) {
  R = R0;
  a = a0;
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  Wstart = Wstart0;
  Efermi = Efermi0;
  beta = beta0;
  EnergyType = EnergyType0;
  form.init(1.,R,a);
  TF.init(A,B,C,D,Wstart,Efermi);

}
//***************************************************************
  /**
   * load parameter where the totla strength is specified not the 
   * energy dependence. THis is useful for OM fits to single energy data
   \param W0 is the strength of the surface imaginary potential
   \param R0 is the radius of the potential in fm
   \param a0 is the diffuseness of the potential in fm
   \param beta0 is the nonlocality length in fm
  */
void surfaceTF::load(double W0, double R0, double a0, double beta0)
{
  R = R0;
  a = a0;
  beta = beta0;
  W = W0;
  form.init(1.,R,a);
  EnergyType = 2;

}
//***********************************************************
  /**
   * for each given energy this needs to be runs before calling the 
   *subsequent functions with give radial dependencies.
   */

void surfaceTF::setEnergy(double Ecm0) {

    Ecm = Ecm0;

    if ( EnergyType == 0 ) {

        if ( Ecm <= Efermi ) W = TF.funct( Ecm ); 
        else W = 0;

        TF.deltaVXAboveBelow( Ecm - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        V = TF.deltaVBelow(); // Subtracted Dispersive Correction
        derV = TF.derDeltaVBelow(); // Derivative

    }
    else if ( EnergyType == 1 ) {

        if ( Ecm >= Efermi ) W = TF.funct( Ecm );
        else W = 0;

        TF.deltaVXAboveBelow( Ecm - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        V = TF.deltaVAbove(); // Subtracted Dispersive Correction
        derV = TF.derDeltaVAbove(); // Derivative
    }
    else {

        W = TF.funct(Ecm);
        V = TF.deltaV(Ecm);
        derV = TF.derDeltaV( Ecm );

    }

}
//***************************************************************
double surfaceTF::imaginaryE( double E ) {

    if ( EnergyType == 0 ) {

        if ( E <= Efermi ) return TF.funct( E ); 
        else return 0;

    }
    else if ( EnergyType == 1 ) {

        if ( E >= Efermi ) return TF.funct( E );
        else return 0;

    }
    else return TF.funct(Ecm);

}

//***************************************************************
double surfaceTF::dispersiveE( double E ) {

    if ( EnergyType == 0 ) {

        TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        return TF.deltaVBelow(); // Subtracted Dispersive Correction
    }
    else if ( EnergyType == 1 ) {

        TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        return TF.deltaVAbove(); // Subtracted Dispersive Correction
    }
    else return TF.deltaV(Ecm);

}
//***************************************************************
double surfaceTF::derDispersiveE( double E ) {

    if ( EnergyType == 0 ) {

        TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        return TF.derDeltaVBelow(); // Derivative
    }
    else if ( EnergyType == 1 ) {

        TF.deltaVXAboveBelow( Ecm - Efermi ); // Dispersive Correction
        TF.findHalfIntegral(); // Needed for Subtracted Dispersive Correction

        return TF.derDeltaVAbove(); // Subtracted Dispersive Correction
    }
    else return TF.derDeltaV(Ecm);
}

//***************************************************************
/**
 * returns the surface imaginary potential at a given radius.
 * the fubction setEnergy(Ecm) must be run before using this function
\param r is radial distance in fm
 */
double surfaceTF::ImaginaryPot(double r)
{
  return W*form.DerWoodSaxon(r);
}
//***************************************************************
  /**
   *returns the real dispersive correction to the surface imaginary
   *potential. The function setEnergy(Ecm) must be run beforehand
   \param r is the radial distance in fm
  */
double surfaceTF::DispersiveCorrection(double r)
{
  return V*form.DerWoodSaxon(r);
}

//***************************************************************
double surfaceTF::DerDispersiveCorrection(double r)
{
  return derV*form.DerWoodSaxon(r);
}
//********************************************
  /**
   * returns the magnitude of the imaginary potential
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceTF::CentralImaginaryPotential(double Ecm)
{
  return TF.funct(Ecm);
}
//*******************************************
  /**
   * returns the magnitude of the dispersive correction
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceTF::CentralDeltaRealPotential(double Ecm)
{
  return TF.deltaV(Ecm);
}
//***************************************
  /**
   *returns the maximum of the imagimary potential as a function of Ecm
   * units are MeV
  */
double surfaceTF::getMaxW()
{
  double Wmax = 0.;
  for (int i=5;i<50;i++)
    {
      double Ecm = (double)i;
      double W = TF.funct(Ecm);
      if (W > Wmax) Wmax = W;
    }
  return Wmax;
}
//******************************************************
double surfaceTF::H(double dr)
{
  if (beta == 0)
    {
      if (dr == 0.) return 1.;
      else return 0.;
    }
  return konst*exp(-pow(dr/beta,2))/pow(beta,3);
}
//******************************************************
complex<double> surfaceTF::U(double r1, double r2)
{
  if (typeN == 0)
    return  -complex<double>(V,W)*
              sqrt(form.DerWoodSaxon(r1)*form.DerWoodSaxon(r2));
  else if (typeN == 1)
     return complex<double>(V,W)*form.DerWoodSaxon((r1+r2)/2.);

  else 
    {
    cout << "typeN = " << typeN << " not known" << endl;
    abort();
    }
}
//******************************************************
// energy derivative of dispersive part
double surfaceTF::dU(double r1, double r2) {

    if (typeN == 0) {

        return  - derV * sqrt(form.DerWoodSaxon(r1)*form.DerWoodSaxon(r2));
    }
    else if (typeN == 1) {

        return derV * form.DerWoodSaxon((r1+r2)/2.);
    }
    else {

        cout << "typeN = " << typeN << " not known" << endl;
        abort();
    }
}
//*******************************************************
complex<double> surfaceTF::potential(double r)
{
  return complex<double>(V,W)*form.DerWoodSaxon(r);
}
//*******************************************************
complex<double> surfaceTF::potential(double r1, double r2, double dr)
{
  return U(r1,r2)*H(dr);
}
//********************************************************
complex<double> surfaceTF::potentialE(double r, double Ecm)
{
  setEnergy(Ecm);
  return potential(r);
}
//*******************************************************
complex<double> surfaceTF::potentialE(double r1, double r2, double dr, 
                                      double Ecm)
{
  setEnergy(Ecm);
  return potential(r1,r2,dr);
}
//*********************************************************
void surfaceTF::setTypeN(int typeN0)
{
  typeN = typeN0;
}
//********************************************************
surfaceTF::surfaceTF()
{
  typeN = 0;
}
//*****************************************************
surfaceTF::surfaceTF(int typeN0)
{
  typeN = typeN0;
}
