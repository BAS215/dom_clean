#include "surfaceGeneral.h"

double const surfaceG::konst  = 1./pow(acos(-1.),3./2.);

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
\param EformType - 0 for absorption that changes smoothly across Efermi, 1 for otherwise.
 */

void surfaceG::load(double R0, double a0, double A0, double B0, double C0,
                     double D0, double Wstart0, double Efermi0,
                     double beta0, int EnergyType0, int EformType0 ) {
  R = R0;
  a = a0;
  A = A0;
  B = B0;
  C = C0;
  D = D0;
  Wstart = Wstart0;
  Efermi = Efermi0;
  betas = beta0;
  EnergyType = EnergyType0; 
  EformType = EformType0;
  form.init(1.,R,a);

  TF.init(A,B,C,D,Wstart,Efermi);

  // form of energy dependence that goes smoothly to zero at 
  // the Fermi energy
  double Ep1 = Wstart;
  int m1 = 4;
  int m2 = 2;
  Eform1.init( A, B, 0.0, Ep1, Efermi, m1, 0, 0., 0. );

  double Ep2 = Ep1 + C;
  double A2 = -A;
  Eform2.init( A2, D, 0.0, Ep2, Efermi, m2, 0, 0., 0. );

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
void surfaceG::load(double W0, double R0, double a0, double beta0)
{
  R = R0;
  a = a0;
  betas = beta0;
  W = W0;
  form.init(1.,R,a);
  EnergyType = 2;

}
//***********************************************************
  /**
   * for each given energy this needs to be runs before calling the 
   *subsequent functions with give radial dependencies.
   */

void surfaceG::setEnergy(double Ecm0) {

    Ecm = Ecm0;

    W = imaginaryE( Ecm );
    V = dispersiveE( Ecm );
    derV = derDispersiveE( Ecm );

//cout<<"Ecm = " <<Ecm <<" E typ and form" <<EnergyType<<","<<EformType <<"\t"<<"Surface real  = "<<V<<"\t"<<"Surface imag =  "<<W<<endl;
}
//***************************************************************
double surfaceG::imaginaryE( double E ) {

    double W0 = 0;
    // Form of energy dependence used in 2007 DOM paper
    if ( EformType == 0 ) {

        if ( EnergyType == 0 ) {

            if ( E <= Efermi ) { 

                W0 = Eform1.ImaginaryPotential( E ) 
                   + Eform2.ImaginaryPotential( E ); 
            }
        }
        else if ( EnergyType == 1 ) {

            if ( E >= Efermi ) { 

                W0 = Eform1.ImaginaryPotential( E )
                   + Eform2.ImaginaryPotential( E );
            }
        }
        else {

            W0 = Eform1.ImaginaryPotential( E ) 
               + Eform2.ImaginaryPotential( E );
        }
    }
    else { // Form of energy dependence used in 2011 DOM paper

        if ( EnergyType == 0 ) {

            if ( E <= Efermi ) W0 = TF.funct( E ); 
        }
        else if ( EnergyType == 1 ) {

            if ( E >= Efermi ) W0 = TF.funct( E );
        }
        else W0 = TF.funct(E);
    }
    return W0;

}

//***************************************************************

double surfaceG::dispersiveE( double E ) {
    double V0;

    // Form of energy dependence used in 2007 DOM paper
    if ( EformType == 0 ) {

        if ( EnergyType == 0 ) {

            // Subtracted Dispersive Correction of imaginary potential
            // below Efermi (even though function says 'Particle')
            V0 = Eform1.DeltaRealParticle( E ) 
               + Eform2.DeltaRealParticle( E ); 
        }
        else if ( EnergyType == 1 ) {

            // Subtracted Dispersive Correction of imaginary potential
            // above Efermi (even though function says 'Hole')
            V0 = Eform1.DeltaRealHole( E ) 
               + Eform2.DeltaRealHole( E ); 
        }
        else {

            // Subtracted Dispersive Correction
            V0 = Eform1.DeltaRealPotential( E ) 
               + Eform2.DeltaRealPotential( E );
        }
    }
    else { // Form of energy dependence used in 2011 DOM paper
        if ( EnergyType == 0 ) {

            TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
            TF.findHalfIntegral(); // Needed for Subtracted Disp. Correction

            V0 = TF.deltaVBelow(); // Subtracted Dispersive Correction
        }
        else if ( EnergyType == 1 ) {

            TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
            TF.findHalfIntegral(); // Needed for Subtracted Disp. Correction

            V0 = TF.deltaVAbove(); // Subtracted Dispersive Correction
        }
        else V0 = TF.deltaV(E);
    }

    return V0;
}
//***************************************************************
double surfaceG::derDispersiveE( double E ) {

    double derV0;
    if ( EformType == 0 ) {

        if ( EnergyType == 0 ) {

            derV0 = Eform1.DerDeltaParticle( E )
                  + Eform2.DerDeltaParticle( E ); // Derivative
        }
        else if ( EnergyType == 1 ) {

            derV0 = Eform1.DerDeltaHole( E )
                  + Eform2.DerDeltaHole( E ); // Derivative
        }
        else {

            derV0 = Eform1.DerDeltaRealPotential( E ) 
                  + Eform2.DerDeltaRealPotential( E );
        }
    }
    else {

        if ( EnergyType == 0 ) {

            TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
            TF.findHalfIntegral(); // Needed for Subtracted Disp. Correction
            derV0 = TF.derDeltaVBelow(); // Derivative
        }
        else if ( EnergyType == 1 ) {

            TF.deltaVXAboveBelow( E - Efermi ); // Dispersive Correction
            TF.findHalfIntegral(); // Needed for Subtracted Disp. Correction
            derV0 = TF.derDeltaVAbove(); // Subtracted Dispersive Correction
        }
        else derV0 = TF.derDeltaV(E);
    }

    return derV0;
}

//***************************************************************
/**
 * returns the surface imaginary potential at a given radius.
 * the fubction setEnergy(Ecm) must be run before using this function
\param r is radial distance in fm
 */
double surfaceG::ImaginaryPot(double r)
{
  return W*form.DerWoodSaxon(r);
}
//***************************************************************
  /**
   *returns the real dispersive correction to the surface imaginary
   *potential. The function setEnergy(Ecm) must be run beforehand
   \param r is the radial distance in fm
  */
double surfaceG::DispersiveCorrection(double r)
{
  return V*form.DerWoodSaxon(r);
}

//***************************************************************
double surfaceG::DerDispersiveCorrection(double r)
{
  return derV*form.DerWoodSaxon(r);
}
//********************************************
  /**
   * returns the magnitude of the imaginary potential
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceG::CentralImaginaryPotential(double Ecm)
{
  return TF.funct(Ecm);
}
//*******************************************
  /**
   * returns the magnitude of the dispersive correction
     \param Ecm in the center-of-mass energy in MeV
   */
double surfaceG::CentralDeltaRealPotential(double Ecm)
{
  return TF.deltaV(Ecm);
}
//***************************************
  /**
   *returns the maximum of the imagimary potential as a function of Ecm
   * units are MeV
  */
double surfaceG::getMaxW()
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
double surfaceG::H(double dr)
{
  if (betas == 0)
    {
      if (dr == 0.) return 1.;
      else return 0.;
    }
  return konst*exp(-pow(dr/betas,2))/pow(betas,3);
}
//******************************************************
complex<double> surfaceG::Un()
{
     return complex<double>(V,W);
}

complex<double> surfaceG::U(double r1, double r2)
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
double surfaceG::dUn() {return derV;}
double surfaceG::dU(double r1, double r2) {

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
complex<double> surfaceG::potential(double r)
{
  return complex<double>(V,W)*form.DerWoodSaxon(r);
}
//*******************************************************
complex<double> surfaceG::potential(double r1, double r2, double dr)
{
  return U(r1,r2)*H(dr);
}
//********************************************************
complex<double> surfaceG::potentialE(double r, double Ecm)
{
  setEnergy(Ecm);
  return potential(r);
}
//*******************************************************
complex<double> surfaceG::potentialE(double r1, double r2, double dr, 
                                      double Ecm)
{
  setEnergy(Ecm);
  return potential(r1,r2,dr);
}
//*********************************************************
void surfaceG::setTypeN(int typeN0)
{
  typeN = typeN0;
}
//********************************************************
surfaceG::surfaceG()
{
  typeN = 0;
}
//*****************************************************
surfaceG::surfaceG(int typeN0)
{
  typeN = typeN0;
}
