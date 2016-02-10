#include "volume.h"

double const volume::konst  = 1./pow(acos(-1.),3./2.);

/**
 *loads all the parameters needed to calculate potentials in the DOM
 *formalism
 */
void volume::load(double Rzero0, double deltaR0, double expR0, 
                 double a0, double A0,
		 double B0, double Ep0, double Efermi0, 
		  int m0, int Asy0, double alpha0, double Ea0, double beta0 ,int Energytype0)
{
  Rzero = Rzero0;
  deltaR = deltaR0;
  expR = expR0;
  a = a0;
  A = A0;
  B = B0;
  Ep = Ep0;
  Efermi = Efermi0;
  m = m0;
  Asy = Asy0;
  alpha = alpha0;
  Ea = Ea0;
  beta = beta0;
  Energytype = Energytype0;
  form.init(1.,Rzero,a);
  volumeF.init(A,B,0.,Ep,Efermi,m,Asy,alpha,Ea);

  if (deltaR > 0.)
    surface.init(A*deltaR/4./a,B,expR,(double)m,Efermi,Ea,alpha,Ep);
  else surface.init(0.,10.,10.,4.,Efermi,Ea,alpha,Ep);
} 
//**********************************************************
//**********************************************************
  /**
   * loads simple imaginary potential of const strength, for use in fitting
   * data from one energy in standard OM
   \param W0 is the strength of the volume imaginary potential in MeV
   \param R0 is the radius in fm
   \param a0 is the diffuseness in fm
   \param beta0 is the nonlocality length in fm
  */
void volume::load(double W0, double R0, double a0, double beta0)
{
  Wvol = W0;
  Vvol = 0.;
  deltaR = 0.;
  Rzero = R0;
  a = a0;
  beta = beta0;
  form.init(1.,Rzero,a);
}
//***********************************************************
  /**
   * in the DOM formalism for each given energy this needs to be run
   * before calling the 
   * subsequent functions with give radial dependies
   \param Ecm0 is the center -f mass energy in MeV
   */
void volume::setEnergy(double Ecm0)
{
  Ecm = Ecm0;
  
  if (Energytype == 0) {

	  if (Ecm <= Efermi) Wvol = volumeF.ImaginaryPotential(Ecm);
	  else Wvol=0;  
	  Vvol = volumeF.DeltaRealParticle(Ecm);
          derVvol = volumeF.DerDeltaParticle( Ecm );
   }
   else if (Energytype == 1) {
	  
          if (Ecm >= Efermi) Wvol = volumeF.ImaginaryPotential(Ecm);
	  else Wvol=0;  
	  Vvol = volumeF.DeltaRealHole(Ecm);
          derVvol = volumeF.DerDeltaHole( Ecm );
   }
   
   else { 
   	
	  
	  Wvol = volumeF.ImaginaryPotential(Ecm);
	  Vvol = volumeF.DeltaRealPotential(Ecm);
          derVvol = volumeF.DerDeltaRealPotential( Ecm );
   }
   
   
  //double R = Rzero;// + deltaR/2.*exp(-pow(B/expR,2)/2.);
  if (deltaR > 0.)
    {
     Wsur = surface.funct(Ecm);

     Vsur = surface.deltaV(Ecm);
     derVsur = surface.derDeltaV( Ecm );

    }

// cout <<"Ecm = " <<Ecm<<" E-type"<<Energytype<<"\t"<<"Wvol = " <<Wvol<< "\t"<< "Vvol = " << Vvol<<endl;
// cout <<"Ecm = " <<Ecm<<"\t"<<"Wvol = " <<Wvol<< "\t"<< "Vvol = " << Vvol<<endl;
}
//***************************************************************
// returns the volume imaginary potential at a given radius
double volume::ImaginaryPot(double r)
{
  double out = Wvol*form.WoodSaxon(r);
  if (deltaR > 0.) out += Wsur*form.DerWoodSaxon(r);
  return out;
}
//***************************************************************
//returns the real dispersive correction to the volume imaginary
//potential
double volume::DispersiveCorrection(double r)
{
  double out = Vvol*form.WoodSaxon(r);
  if (deltaR > 0.) out  +=  Vsur*form.DerWoodSaxon(r);  
  return out;
}
//******************************************************
double volume::DerDispersiveCorrection(double r)
{
  double out = derVvol*form.WoodSaxon(r);
  if (deltaR > 0.) out  +=  derVsur*form.DerWoodSaxon(r);  
  return out;
}
//*********************************************************
void volume::setTypeN(int typeN0)
{
  typeN = typeN0;
}
//********************************************************
volume::volume()
{
  typeN = 0;
}
//*****************************************************
volume::volume(int typeN0)
{
  typeN = typeN0;
}
//******************************************************
double volume::H(double dr)
{
  if (beta == 0.)
    {
      if (dr == 0.) return 1.;
      else return 0.;
    }
  return konst*exp(-pow(dr/beta,2))/pow(beta,3);
}
//******************************************************


complex<double> volume::Un()
{
  complex<double> out;
  out = complex<double>(Vvol,Wvol);
  return out;
}
double volume::dUn() {return derVvol;}
complex<double> volume::U(double r1, double r2)
{
  complex<double> out;
  if (typeN == 0)
    {

      if (deltaR == 0.)out = -complex<double>(Vvol,Wvol)*
	sqrt(form.WoodSaxon(r1)*form.WoodSaxon(r2));
      else if (deltaR > 0.)
	{
	  double vfr = form.WoodSaxon(r1);
          double sfr = form.DerWoodSaxon(r1);
	  double rf1= vfr + Vsur/Vvol*sfr;
          double if1= vfr + Wsur/Wvol*sfr;
           
	  vfr = form.WoodSaxon(r2);
          sfr = form.DerWoodSaxon(r2);
	  double rf2= vfr + Vsur/Vvol*sfr;
          double if2= vfr + Wsur/Wvol*sfr;

          out = -complex<double>(Vvol*sqrt(rf1*rf2),Wvol*sqrt(if1*if2));

	}

    }
  else if (typeN == 1)
    {
      double rr = (r1+r2)/2.;
       out = complex<double>(Vvol,Wvol) * form.WoodSaxon(rr);
	
      if (deltaR > 0.) out += complex<double>(Vsur,Wsur)*
	form.DerWoodSaxon(rr);
    } 
  else 
    {
    cout << "typeN = " << typeN << " not known" << endl;
    abort();
    }


  return out;
}
//******************************************************
// energy derivative 
double volume::dU(double r1, double r2) {

    double out;
    if (typeN == 0) {

        if (deltaR == 0.) {
            out = - derVvol * sqrt(form.WoodSaxon(r1)*form.WoodSaxon(r2));
        }
        else if (deltaR > 0.) {

	        double vfr = form.WoodSaxon(r1);
            double sfr = form.DerWoodSaxon(r1);
	        double rf1= vfr + Vsur/Vvol*sfr;
            double if1= vfr + Wsur/Wvol*sfr;
           
	        vfr = form.WoodSaxon(r2);
            sfr = form.DerWoodSaxon(r2);
	        double rf2= vfr + Vsur/Vvol*sfr;
            double if2= vfr + Wsur/Wvol*sfr;

            out = -derVvol*sqrt(rf1*rf2);

	    }

    }
    else if (typeN == 1) {

        double rr = (r1+r2)/2.;
        out = derVvol*form.WoodSaxon(rr);
        if (deltaR > 0.) out += derVsur * form.DerWoodSaxon(rr);
    } 
    else {

        cout << "typeN = " << typeN << " not known" << endl;
        abort();
    }


    return out;
}
//******************************************************
complex<double> volume::potential(double r)
{
  double ff = form.WoodSaxon(r);
  complex<double> out = complex<double>(Vvol*ff,Wvol*ff);
  if (deltaR > 0.) 
    {
      ff = form.DerWoodSaxon(r);
      out += complex<double>(Vsur*ff,Wsur*ff);
    }
  return out;
}
//**********************************************************
complex<double> volume::potential(double r1, double r2, double dr)
{
  return U(r1,r2)*H(dr);
}
//*************************************************************
complex<double> volume::potentialE(double r, double Ecm)
{
  setEnergy(Ecm);
  return potential(r);
}
//**************************************************************
complex<double> volume::potentialE(double r1, double r2, double dr, double Ecm)
{
  setEnergy(Ecm);
  return potential(r1,r2,dr);
}
//****************************************************
/**
 * sets the equivalent local potential as per 
 * Perey and Buck
\param WvolE0 is depth of equivalent local volume potential
*/
void volume::setEquivLocalDepth(double WvolE0)
{
  WvolE = WvolE0;
}
//****************************************************
  /**
   * return the equivalent local potential as
   * per Perey and Buck
    \param r - radius in fm
   */
double volume::potentialEquivLocal(double r)
{
  return WvolE*form.WoodSaxon(r);
}
