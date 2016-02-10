#include "spinOrbit.h"

/**
 * loads the DOM spin orbit parameters
 */
void spinOrbit::load(double R0, double a0, 
     double Vzero0, double Efermi0,double AW0, double BW0)
{
  R = R0;
  a = a0;
  Vzero = Vzero0;
  Efermi = Efermi0;
  AW = AW0;
  BW = BW0;
  Form.init(AW,BW,0.,0.,Efermi,4,0,0.,0.);
  
}
//**********************************************************
/**
 * loads the strength of the real and imaginary spinorbit
 * this version of load is udeful for OM fits to single-energy data
 \param V0 is real strength of spin orbit potential
 \param W0 is imaginary stregth of spin orbit potential
 \param R0 is radius of spin-orbit potential in fm
 \param a0 is the diffuseness of spin-orbit potential in fm
*/
void spinOrbit::load(double V0, double W0, double R0, double a0)
{
  R = R0;
  a = a0;
  V = V0;
  W = W0;
  Real.init(V,R,a);
  Disp.init(V,R,a);
  DerDisp.init( V, R, a );
  Imag.init(W,R,a);

}

//*******************************************************************
void spinOrbit::setEnergy(double Ecm0)
{
  Ecm = Ecm0;

  double V_disp= Form.DeltaRealPotential( Ecm );
  double der_V_disp = Form.DerDeltaRealPotential( Ecm );

 //V_disp = 0 for MSU business
  V  = Vzero + V_disp;

  Real.init(V,R,a);
  Disp.init( V_disp, R, a ); //Disperisve correction sjw 02/21/2012
  DerDisp.init( der_V_disp, R, a ); // Energy derivative sjw

  W = Form.ImaginaryPotential(Ecm);
  Imag.init(W,R,a);
  //cout<<"Ecm =" <<Ecm <<"\t"<< "SO real =  " <<V << "\t SO imag" << W <<"\t V_disp =  "<<V_disp<<"   Vzero= "<< Vzero<<endl;

}
//***************************************************************
double spinOrbit::RealPotential(double r)
{
  return Real.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}

//***************************************************************
double spinOrbit::ImaginaryPotential(double r)
{
  return Imag.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}

//***************************************************************
double spinOrbit::DispersiveCorrection(double r)
{
  return Disp.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}

//***************************************************************
double spinOrbit::DerDispersiveCorrection(double r)
{
  return DerDisp.DerWoodSaxon(r)/r/2.*LdotSigma/a;
}

//****************************************************************
complex<double> spinOrbit::potential(double r)
{
  return complex<double>(RealPotential(r),ImaginaryPotential(r));
}
//****************************************************************
complex<double> spinOrbit::potentialE(double r, double Ecm)
{
  setEnergy(Ecm);
  return complex<double>(RealPotential(r),ImaginaryPotential(r));
}
//*****************************************************************
double spinOrbit::setAM(int l, double j)
{
  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  return LdotSigma;
}
