#include "boundRspace.h"

double const boundRspace::pi =  acos(-1.);
double const boundRspace::kconstant =  .048192;
double const boundRspace::e2 = 1.44; // Coulomb constant

void boundRspace::init( double rmax0, int hamiltonian_pts0, double Efermi0, 
                        double ph_gap0, int Lmax0, double Z0, double Zp0, 
                        double A0) {

  rmax = rmax0;
  hamiltonian_pts = hamiltonian_pts0;
  Efermi = Efermi0;
  ph_gap = ph_gap0;
  Lmax = Lmax0;
  Z = Z0;
  A = A0;
  Zp = Zp0;
 // mu = 1.;
  mu = A/(A-1.);
}
//*******************************************************
void boundRspace::newEnergy(double Ecm0)
{
  Ecm = Ecm0;



  Kwave2 =  kconstant*mu*Ecm;

  Kwave = sqrt(abs(Kwave2));

  muhbar = kconstant*mu;

  gamma = fabs((Z-1)*Zp*e2*Kwave/Ecm/2.); //Sommerfeld parameter

  Pot->setEnergy(Ecm);

}
//*********************************************************
boundRspace::boundRspace( double rmax0, int hamiltonian_pts0, double Efermi0,
                          double ph_gap0, int Lmax0, double Z0, double Zp0,
                          double A0, pot*Pot0 ) {

  init( rmax0, hamiltonian_pts0, Efermi0, ph_gap0, Lmax0, Z0, Zp0, A0 );
  Pot = Pot0;

  // radial grid for Numerov code
  rStart = .05;
  rStop = 12.;
 // Nnumerov = 200; 
  Nnumerov = 200; 

  rdelt = rmax / hamiltonian_pts;

  mWave = Nnumerov + 1;
  nWave = mWave + 120;  

  WaveArray = new double [nWave];

 // initialize the coupled differential equation solver
  initNumerovR(rStart,rStop,Nnumerov);

  // misc
  max_n = 10; // maximum number of principal quantum numbers

  //set up array of Log derivatives.
  Nl = 7;
  LogDerMin = -123;
  LogDerMax = 5;
  NLogDer = LogDerMax - LogDerMin + 1;
  
  LogDerArray = new double *[Nl+1];
  for (int i=0;i<=Nl;i++) LogDerArray[i] = new double [NLogDer];


 // cout << " Making LogDer array" << endl;
  for (int l = 0;l<=Nl;l++)
    {

     for (int i=LogDerMin;i<=LogDerMax;i++)
      {
        Ecm = (double) i;
        Kwave2 = kconstant*mu*Ecm;
        muhbar = kconstant*mu;
        Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1
        gamma = mu*(Z-1)*Zp/Kwave/28.820;  // Sommerfeld parameter

        LogDerArray[l][i-LogDerMin] = exteriorLogDer(l);
        
      }
    }
//  cout << " Finished LogDer array" << endl;


}
//*********************************************************************
boundRspace::~boundRspace()
{
  delete [] WaveArray;

  for (int i=0;i<=Nl;i++) delete [] LogDerArray[i];
  delete [] LogDerArray;
}
//****************************************************************************
  /**
   * scans energy to find a bound state
   \param Elower - starting energy in Mev
   \param Eupper - finishing energy in MeV
   \param j - total spin a level
   \param l - orbital AM of level
    \param Ifine - =1 fine search otherwise less accurate
   */

int boundRspace::searchNonLoc(double& Elower,double Eupper, double j, int l, int Ifine)
{

  int iout;
  phase = 0;  // equivalent local potential from Perey and Buck 
  iout = searchLoc(Elower,Eupper,j,l,Ifine); // find bound state with this potential
  Elower0 = Elower;

  //correct local wavefunction for nonlocality

  if (Elower < 0.)
    {
    //correct the local wavefunction for nonlocality
    //the correction at r=0, must be unity, so I scale the usual correction
    //to make this so.
    Pot->HartreeFock.setEquivLocalDepth(1.- sqrt(mstar));
    for (int i=0;i<=Nnumerov+1;i++) y[i] *=
      (1. +  Pot->HartreeFock.potentialEquivLocal(x[i]))/sqrt(mstar); 
    }

  //if (l == 3 && j == 3.5) cout << "start " << 0 << " " << l << " " << j << " " << Elower << " " << iout << endl;
  if (iout == 0) return iout;


  phase = 1; // use TELP equiv local potential
  generateVloc = false; // do not create a new Equiv Local poential
                                 //at each energy  
   // for (int itry = 1;itry<4;itry++)
    for (int itry = 1;itry<6;itry++)
      {
	if ( l == 0 && itry == 1)Elower -= 15.;
        else if ( l == 1 && j == 0.5 && itry == 1 ) Elower -=10.;
	else if ( l == 1 && j == 1.5 && itry == 1) Elower -= 10.;
	else if (itry == 1)Elower -= 7.;
	else Elower -= 5.;

	generateLocEquiv(l); //generate new equivalent local pot( TELP)


        iout = searchLoc(Elower,Eupper,j,l,Ifine); // find bound state with this potential
	//if (l == 3 && j == 3.5)cout << "iter= " << itry << " " <<  l << " " << j << " " << Elower << " " << iout << endl;
	if (iout == 0) return iout; 

      }
    return iout;
}
//***************************************
//***************************************************************
  /** determines the bound state and quasibound states energies
   * for local potential
    * if one wants to get the wavefunction, you need to get the energy
    *  much more accurately than this, use Ifine=1
    * search is made in the energy range Elower to Eupper.
    * found level is returned as Elower.
    *function returns 1 is level is found, else 0.
    \param Elower is the lower bound for energy range of search  [MeV]
    \param Eupper is the upper bound for energy range of search  [MeV]
    \param j total angular momentum of the single-particle state
    \param l is the orbital angular momentum of the single-particle state
    \param Ifine fine search of energy - need if you want wavefunctions
  */
int boundRspace::searchLoc(double& Elower,double Eupper, double j, int l, int Ifine)
{

  double const deltaE = 1;
  double Ecm  = Elower;

  double Wave;
  double WaveOld = 0.;
  double WaveAbove = 0.;
  double EcmAbove = 0.;
  int tries = 0;
  int last = 0;



  for (;;)
    {
      // find all energy dependent OM terms and initialize OM
      newEnergy(Ecm); 



       // integrate wavefunction and find phaseshifts

       integrateWaveFunction(j,l);
       Wave = LogDerDifference(j,l);
       //if ( l == 3 && j == 3.5) cout << Ecm << " " << Wave << endl;
       if (last) 
	 {
           if (Ifine)
	     {

	       if (Wave*WaveOld < 0.) 
		 {
		   Elower = Ecm - deltaE/2.;
                   Eupper = Ecm;
		 }
               else 		 {
                   Elower = Ecm;
                   Eupper = EcmAbove;
		 }
               return searchFine(Elower,Eupper,j,l);

	     }
	   else  // last iteration, use Ridder method for an interpolation
	 //see Numerical Recipes in Fortran77 Page351
	     {
               double sign = 1.;
               if (WaveOld < WaveAbove) sign = -1.;
               Elower = Ecm + deltaE/2.*sign*Wave/
	        sqrt(pow(Wave,2)-WaveOld*WaveAbove);
               return 1;
	     }
	 }

       if (tries > 0)
	 {
	   if (Wave*WaveOld < 0.)
	     {
               EcmAbove = Ecm;
               Ecm -= deltaE/2.;
               WaveAbove = Wave;
	       last = 1; // flag the last iteration
	       continue;
	     }
	 }

       tries++;
       WaveOld = Wave;
       Ecm += deltaE;
       if (Ecm > Eupper) break;
 
    }
  return 0;
}

//*********************************************************************
  /**
   * After integration the wavefunction to rStop, this function returns 
   *Wavefunction magnitude and its derivative
   \param l0 is the orbital angular momentum of the nucleon
   \param j0 is the total angular momentum of the nucleon
   */
void boundRspace:: integrateWaveFunction( double j, int l)
{
 

  Pot->SpinOrbit.setAM(l,j);

  double y0 = pow(rStart,l+1);
  double dydr0 = (double)(l+1)*pow(rStart,l);
  
  if (phase == 0)
    {
      nonLocalArray0(l); //use perey and Buck equiv local potential(1st iteration)
     solveNumerovR(y0,dydr0);
    }
  else if (phase == 1)
    {
      nonLocalArray1(false,l); //use TELP equiv local potential(subsequent inerations)
     solveNumerovR(y0,dydr0);
    }
  else 
    integrateWaveFunction0(y0,dydr0,l);



  y_end =   y[Nnumerov];
  dydr_end = (y[Nnumerov+1]-y[Nnumerov-1])/2./dx;

}
//****************************************************************
double boundRspace::LogDerDifference( double j, int l)
{
  //this would be fine if you didn't need the wavefunction
  //valarray<double> Waveff = IntegrateBound(j,l);



  //log derivative of interior wave function
  //double LogDerInterior = Waveff[1]/Waveff[0];

  //interpolate the exterior LogDerative
  double en = floor(Ecm);
  double delta = Ecm-en;
  int i = (int)en- LogDerMin;
  if (i < 0 || i > NLogDer) cout << " outside of LogDerArray, Ecm = " 
				<< Ecm << " i= " << i << endl;
  double LogDerExterior = LogDerArray[l][i];

  if (delta != 0.)
    {
      if (Zp == 0 && fabs(Ecm)< 1.) 
	LogDerExterior = exteriorLogDer(l);

      else LogDerExterior += delta*(LogDerArray[l][i+1]-LogDerExterior);
    }

  // if (l == 4) cout << "ext = " << LogDerExterior << endl;

  //tied the difference in interior and exterior logderivatives, 
  //but is not very very useful
  //as it only changes very close to a bound state andf hence they easily
  //missed.
  //return LogDerInterior - LogDerExterior;

  //instead
  //normalise interior and exterior wavefunctions to have the same slope,
  //then take the different in the magnitude of these wavefunctions - 
  // this crosses zero at a eigenfunction. one can easily find the zeros.

  return y_end - dydr_end/LogDerExterior;
  
}
//************************************************************************
  /**
   * finds the eigenstates accuracy enough so that the Wavefunction
   *can be determined to give physical properties of the state.
   *Ridder method for an interpolation,
   *see Numerical Recipes in Fortran77 Page351
    \param Elower is the lower bound for energy range of search  [MeV]
    \param Eupper is the upper bound for energy range of search  [MeV]
    \param j total angular momentum of the single-particle state
    \param l is the orbital angular momentum of the single-particle state
   */
int boundRspace::searchFine(double& Elower,double Eupper, double j, int l)
{
  int counter = 0;
  if (Elower >= Eupper) 
  cout << "Hey buddies, what gives? Elower >= Eupper" << endl;

  double Ecm;

  // start at lower limit, find difference in slope
  // of integrated wavefunction and matched whittaker function
  Ecm = Elower;
  newEnergy(Ecm);
  integrateWaveFunction(j,l);
  double ylower = LogDerDifference(j,l);

  //------------------------------------------------------
  // next the upper limit
  Ecm = Eupper;
  newEnergy(Ecm);
  integrateWaveFunction(j,l);
  double yupper = LogDerDifference(j,l);




  if (yupper*ylower > 0.)
    { 
    cout << "hey what gives, zero is not enclosed by limits" << endl;
    cout << ylower << " " << yupper << endl;
    cout << Elower << " " << Eupper << endl;
    }

  for (;;)
    {

      //cout << Elower << " " << ylower << " " << 
      // Eupper << " " << yupper << endl;

      double Emiddle  = (Elower+Eupper)/2.;
      double Ecm  = Emiddle;
 
      //return if precision is too small
      if (Ecm == Elower) return 1;
      if (Ecm == Eupper) return 1;
      newEnergy(Ecm);
       integrateWaveFunction(j,l);
      double ymiddle = LogDerDifference(j,l);
      //if ( l == 3 && j == 3.5)cout << "mid = " << Ecm << " " << ymiddle << endl;

       double sign = -1.;
       if (ylower > yupper) sign = 1.;
       
       Ecm = Emiddle + (Emiddle-Elower)*sign*ymiddle/
	 sqrt(pow(ymiddle,2)-ylower*yupper);

       newEnergy(Ecm);
       integrateWaveFunction(j,l);
       double y = LogDerDifference(j,l); 
       //if (l == 3 && j == 3.5)cout << " new " << Ecm << " " << y << endl;

       if (abs(y) < .00001) 
	 {
	   Elower = Ecm;
           return 1;
	 }

       // find new bounding energies
       if (Ecm < Emiddle)
	 {
           if (y*ylower < 0.)
	     {
	       Eupper = Ecm;
               yupper = y;
	     }
	   else 
	     {
               Elower = Ecm;
               Eupper = Emiddle;
               ylower = y;
               yupper = ymiddle;
	     }
	 }
       else 
	 {
           if (y*yupper < 0.)
	     {
               Elower = Ecm;
               ylower = y;
	     }
	   else
	     {
               Elower = Emiddle;
               Eupper = Ecm;
               ylower = ymiddle;
               yupper = y;
	     }
	 }
       if (counter > 11) return 0;
       counter++;

    }
  return 0;
}

//*********************************************************************
void boundRspace::integrateWaveFunction0(double y0, double dydr0, int l)
  {
    nonLocalArray0(l);
    solveNumerovR(y0,dydr0);

    if (Pot->beta_max == 0.) return;


  //  for (int itry = 1;itry<8;itry++)
    for (int itry = 1;itry<12;itry++)
      {
	nonLocalArray1(true,l);
        solveNumerovR(y0,dydr0);

      }
  }

//**********************************************************



void boundRspace::generateLocEquiv(int l)
{
  for (int i=0;i<=Nnumerov+1;i++)
    {
	
      //get F factor
      double F = 0;
      if (i > 0)
        {
          double slope;

          //estimate slope of wavefunction   
          if (i < Nnumerov+1) slope = (y[i+1]- y[i-1])/2./dx;
	  else slope = (3.*y[i] + y[i-2] - 4.*y[i-1])/2./dx;


          double stuff= y[i]/slope;
      
          F = exp(-100.*pow(stuff,2));
          stuff = pow(x[i],l+1)/y[i] -1.;

          F *= 1. - exp(-100.*pow(stuff,2));
        }

       // now do the convolution
      double v = 0.;
      double r1 = x[i];
      for (int j=0;j<=Nnumerov+1;j++)
        {
	  double r2 = x[j];
	  if (r1 - r2 > 3.*Pot->beta_max) continue;
	  if (r2 - r1 > 3.*Pot->beta_max) break;
	  // double ddr = abs(r1-r2)/3.; 

        v += dx * y[j] * real( Pot->nonlocalPart( r1, r2 ) );
/*
          if (Pot->HartreeFock.beta0 > 0. && ddr < Pot->HartreeFock.beta0)
	    v += dx*Pot->angleIntegration(r1,r2,Pot->HartreeFock.beta0,l)*y[j]*
	      Pot->HartreeFock.U0(r1,r2);

	   if (Pot->HartreeFock.beta1 > 0. && ddr < Pot->HartreeFock.beta1 
               && Pot->HartreeFock.A != 0.)
	     v += dx*Pot->angleIntegration(r1,r2,Pot->HartreeFock.beta1,l)*y[j]*
	       Pot->HartreeFock.U1(r1,r2);

	   if (Pot->Volume.beta > 0. && ddr < Pot->Volume.beta)
	     {
	       double two =  (Pot->Volume.U(r1,r2)).real();
	       double three = (Pot->Surface.U(r1,r2)).real();
	       double fr = two + three;
	       v += dx*Pot->angleIntegration(r1,r2,Pot->Volume.beta,l)*y[j]*fr;
	     }
*/
	 }


       source[i] = muhbar*F*v;
       vloc[i] = muhbar*(1.-F)/y[i]*v;
    }
}
//******************************************************
void boundRspace::nonLocalArray1(bool generateVloc, int l)
{

  //source term and local potentials
  if (generateVloc) generateLocEquiv( l);

  for (int i=0;i<=Nnumerov+1;i++)
    {

      // now add in the remainding local parts of the potential
      //start with coulomb
      double fr = Pot->coulomb(x[i]);


      // next the spin-orbit potential
      fr += Pot->SpinOrbit.RealPotential(x[i]);


      fr += Pot->HartreeFock.potentialSurface(x[i]);

      if (Pot->HartreeFock.beta0 == 0.)
      fr += Pot->HartreeFock.potential0(x[i]);

      if (Pot->HartreeFock.beta1 == 0. && Pot->HartreeFock.A != 0.)
      fr += Pot->HartreeFock.potential1(x[i]);

      
      if (Pot->VolumeAbove.beta == 0.)
	{
         double two = Pot->VolumeAbove.DispersiveCorrection(x[i]);
	 fr += two ;
	}

      
      if (Pot->SurfaceAbove.betas == 0.)
	{
         double three = Pot->SurfaceAbove.DispersiveCorrection(x[i]);
	 fr += three;
	}

      if (Pot->VolumeBelow.beta == 0.)
	{
         double four = Pot->VolumeBelow.DispersiveCorrection(x[i]);
	 fr += four;
	}
      

      if (Pot->SurfaceBelow.betas == 0.)
	{
         double five = Pot->SurfaceBelow.DispersiveCorrection(x[i]);
	 fr +=  five;
	}
       
      fr *= muhbar;
      if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);
      fr = -(Kwave2 - fr); 


      f[i] = vloc[i] + fr;
 
    }
}
//*********************************************************
  /**
   * sets the source and F arrays for the the first
   * iteration of the wavefunction for the nonlocal calculation
   */
void boundRspace::nonLocalArray0(int l)
{
  double coulDis = 1.73*(Z-1)*Zp/Pot->Rc;
  //double V = - Pot->HartreeFock.Vvol - Pot->VolumeAbove.Vvol;

  double V = - Pot->HartreeFock.Vvol - Pot->VolumeAbove.Vvol- Pot->HartreeFock.V_wine;


  double k2;
 // for (int i=0;i<16000;i++)
  for (int i=0;i<100000;i++)
    {
     k2 = (Ecm - V - coulDis)*kconstant*mu;
     double Vnew = -Pot->HartreeFock.Vvol/(1.+Pot->HartreeFock.A)
        *exp(-k2*pow(Pot->HartreeFock.beta0,2)/4.)
        - Pot->HartreeFock.A*Pot->HartreeFock.Vvol/(1.+Pot->HartreeFock.A)*
           exp(-k2*pow(Pot->HartreeFock.beta1,2)/4.)
       - Pot->VolumeAbove.Vvol*exp(-k2*pow(Pot->VolumeAbove.beta,2)/4.)
       - Pot->VolumeBelow.Vvol*exp(-k2*pow(Pot->VolumeBelow.beta,2)/4.);
       - Pot->HartreeFock.V_wine*exp(-k2*pow(Pot->HartreeFock.R_wine,2)/4.);

     V = (V+Vnew)/2.; //this step helps damp an oscillation with each step
                     // of the iteration


    }
  //effective mass relative to nucleon mass
  mstar = 
    20.9205*k2/(Ecm + Pot->HartreeFock.Vvol + Pot->VolumeAbove.Vvol 
                + Pot->VolumeBelow.Vvol - coulDis + Pot->HartreeFock.V_wine);

   V*= 1.1;
  Pot->HartreeFock.setEquivLocalDepth(-V);

  double factorI1 = exp(-k2*pow(Pot->VolumeAbove.beta,2)/4.);
  double factorI2 = exp(-k2*pow(Pot->VolumeBelow.beta,2)/4.);

    for (int i=0;i<N+2;i++)
      {
	double fr = Pot->HartreeFock.potentialEquivLocal(x[i]);

	double two = 0.;//Pot->Volume.DispersiveCorrection(x[i])*factorI;
       double three = Pot->SurfaceAbove.DispersiveCorrection(x[i])*factorI1;
       double four = Pot->SurfaceBelow.DispersiveCorrection( x[i] ) * factorI2;
       fr += two + three + four;
       
       fr += Pot->HartreeFock.potentialSurface(x[i]);

       fr += Pot->coulomb(x[i]);
       fr += Pot->SpinOrbit.RealPotential(x[i]);


       fr *= muhbar;

       if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);
       fr = -(Kwave2 - fr);


       f[i] = fr;

       source[i] = 0.;
      }
}



//*********************************************************************

  /**
   * calculates the log derivative of the exterior wavefunction
   \param l is the orbital angular momentum of the nucleon
  */
double boundRspace::exteriorLogDer(int l)
{
  if (Zp == 0.)
    {
      sphericalB sph;
      if (Ecm < 0.) return sph.LogDer_k(l,rStop*Kwave)*Kwave;
      else if (Ecm >0.) return sph.LogDer_y(l,rStop*Kwave)*Kwave;
      else return -1.e32;
    }
  else
    {

      if (Ecm <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  double outwave,DoutwaveDr;
	  if (Ecm < 0.)
	    {
	      outwave = whittaker.getWaveFunction(gamma,l,rStop*Kwave);
	      DoutwaveDr = whittaker.derivative*Kwave;
	    }
	  else 
	    {
	      if (Zp  == 0.) return 1.e32;
	      outwave = whittaker.ZeroEnergy( 2.*mu*(Z-1)*Zp/28.820*rStop,l);
	      DoutwaveDr = whittaker.derivative*2.*mu*(Z-1)*Zp/28.820;
	    }
	  return DoutwaveDr/outwave;
	}
      else  // from irregular Coulomb wave function for quasi-bound states
	{
	  coul Coulomb;
	  Coulomb.init(l,gamma,rStop*Kwave);
	  return Coulomb.dG/Coulomb.G*Kwave;
	}
    }
}
//****************************************************
  /**
   *determines the number of nodes in the wavefunction
   */
int boundRspace::nodes()
{
  int nodes = 0;
  for (int i = 1;i<mWave-2;i++)
    {
 
      //cout << i << " " <<y[i] << endl;
      if (y[i]*y[i-1] <= 0.) nodes++;

    }
  return nodes;
}
//**************************************************
  /**
   * normalise the radial wavefunction to unity
   */
void boundRspace::normalizeWF()
{

  double sum = 0.;
  for (int i=0;i<nWave;i++)
    {
      sum += pow(WaveArray[i],2);
    }
  sum *= dx;
  sum = sqrt(sum);
  for (int i=0;i<nWave;i++)
    {
      WaveArray[i] /= sum;
    } 
}
//************************************************
  /**
   * adds the exterior part of the wavefunction to the array WaveArray
   * containing the interior part.
   \param l is the orbital angular momentum of the nucleon
  */
void boundRspace::exteriorWaveFunct(int l)
{

  for (int i=0;i<mWave;i++)
    {
      WaveArray[i] = y[i];
    }

  double ANC=0.;

  if (Zp == 0)
    {
      sphericalB sph;
      for (int i=0;i<=nWave-mWave;i++)
	{
	 double r =rStop + (double)i*dx;
	 double outwave;
         if (Ecm <= 0.) outwave = sph.k(l,r*Kwave);
	 else outwave = sph.y(l,r*Kwave);
	 if (i==0) ANC = WaveArray[mWave-1]/outwave;
	 else WaveArray[mWave+i-1] = outwave*ANC;
	}
    }
  else
    {
      if (Ecm <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*dx;
	      double outwave;
	      if (Ecm < 0.) outwave = 
              whittaker.getWaveFunction(gamma,l,r*Kwave);

	      else outwave = whittaker.ZeroEnergy( 2.*mu*(Z-1)*Zp/28.820*r,l);
	      if( i==0) ANC = WaveArray[mWave-1]/outwave;
	      else WaveArray[mWave+i-1] = outwave*ANC;
	    }
	}
      else  // from irregular Coulomb wave function for quasi-bound states
	{
	  coul Coulomb;
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*dx;
	      Coulomb.init(l,gamma,r*Kwave);
	      if( i==0) ANC = WaveArray[mWave-1]/Coulomb.G;
	      else WaveArray[mWave+i-1] = Coulomb.G*ANC;
	    }

	}
    }
}

// *************************************************************
// The following set of functions allow for the calculation of
// bound states by diagonalization of the Hamiltonian.
//
// *************************************************************

// returns the real part of the Hamiltonian in coordinate space.
// This is needed when solving the Schroedinger-like equation 
// for the overlap functions.
//
std::vector< double >
boundRspace::make_rmesh_point( ) {

    std::vector< double > rmesh;
    for ( int i = 0; i < hamiltonian_pts; ++i ) {
       
       // use the midpoint rule for integration
//      if (i==0) {rmesh.push_back( ( i + 0.5 ) * rdelt) ;} 
    //    rmesh.push_back( ( i + 0.5 ) * rdelt * 1.02564 ); 
      rmesh.push_back( ( i + 0.5 ) * rdelt * (39./40.) ); 
    //   rmesh.push_back( ( i + 0.5 ) * rdelt ); 
    }

    return rmesh;
}
std::vector< double >
boundRspace::make_rmesh( ) {

    std::vector< double > rmesh;
    for ( int i = 0; i < hamiltonian_pts; ++i ) {
       
       // use the midpoint rule for integration
       rmesh.push_back( ( i + 0.5 ) * rdelt ); 
    }

    return rmesh;
}

matrix_t 
boundRspace::re_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J ) {

    Pot->setEnergy( Ecm ); 
    Pot->setAM( L, J );

    double rdelt = rmesh.at(1) - rmesh.at( 0 );

    // initialize matrix
    matrix_t ham( rmesh.size(), rmesh.size() );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    for ( unsigned int j = 0; j <= i; ++j ) {
        
        ham( i, j ) = 0.0;
        ham( j, i ) = ham( i, j );
    }
    }

    // \hbar^2 / ( 2 * mu ), mu is reduced mass
    double fac = - 1 / ( Pot->kconstant * mu );

    // Construct Tridiagonal matrix for kinetic energy term
    // and add in diagonal parts of potential
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        // Diagonal Part of Kinetic Energy
        double KE_diagonal = -2 * fac / std::pow( rdelt, 2 )
                          - fac * L * ( L + 1 ) / std::pow( rmesh[i], 2 );

        // Off-diagonal part of Kinetic energy
        double KE_off_diagonal = fac / std::pow( rdelt, 2 );

        // Diagonal Parts of Potential Energy
        double PE_diagonal = real ( Pot->localPart( rmesh[i] ) );

        ham( i, i ) = PE_diagonal + KE_diagonal;
        if ( i != 0 ) ham( i, i - 1 ) = KE_off_diagonal;
        if ( i != rmesh.size() - 1 ) ham( i, i + 1 ) = KE_off_diagonal;

        // Boundary Condition
        if ( i == 0 ) ham( 0, 0 ) -= KE_off_diagonal * std::pow( -1.0, L );

        // Now add in nonlocal components
        for ( unsigned int j = 0; j <= i; ++j ) {

            // nonlocalPart is already multiplied by r, r'
            ham( i, j ) += real( Pot->nonlocalPart( rmesh[i], rmesh[j] ) ) * rdelt;

            if ( j != i ) ham( j, i ) = ham( i, j );

        } // end loop over j

    } // end loop over i

    return ham;
}

// returns the complex part of the Hamiltonian. This is needed when
// solving the Dyson equation to get the propagator 
cmatrix_t 
boundRspace::c_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J) {

    Pot->setEnergy( Ecm ); 
    Pot->setAM( L, J );

    double rdelt = rmesh.at(1) - rmesh.at( 0 );

    // initialize matrix
    cmatrix_t ham( rmesh.size(), rmesh.size() );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    for ( unsigned int j = 0; j <= i; ++j ) {
        
        ham( i, j ) = complex<double>( 0.0, 0.0 );
        ham( j, i ) = ham( i, j );
    }
    }

    // \hbar^2 / ( 2 * mu ), mu is reduced mass
    double fac = - 1 / ( kconstant * mu );

    // Construct Tridiagonal matrix for kinetic energy term
    // and add in diagonal parts of potential
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        // Diagonal Part of Kinetic Energy
        double KE_diagonal = -2 * fac / std::pow( rdelt, 2 )
                          - fac * L * ( L + 1 ) / std::pow( rmesh[i], 2 );

        // Off-diagonal part of Kinetic energy
        double KE_off_diagonal = fac / std::pow( rdelt, 2 );

        // Diagonal Parts of Potential Energy
        complex<double> PE_diagonal = Pot->localPart( rmesh[i] );

        ham( i, i ) = PE_diagonal + KE_diagonal;
        if ( i != 0 ) ham( i, i - 1 ) = KE_off_diagonal;
        if ( i != rmesh.size() - 1 ) ham( i, i + 1 ) = KE_off_diagonal;

        // Boundary Condition
        if ( i == 0 ) ham( 0, 0 ) -= KE_off_diagonal * std::pow( -1.0, L );

        // Now add in nonlocal components
        for ( unsigned int j = 0; j <= i; ++j ) {

            // nonlocalPart is already multiplied by r, r'
            ham( i, j ) += Pot->nonlocalPart( rmesh[i], rmesh[j] ) * rdelt;

            if ( j != i ) ham( j, i ) = ham( i, j );

        } // end loop over j
    } // end loop over i

    return ham;
}


// This actually calculates r * r' * rdelt * G( r, r' )
cmatrix_t 
boundRspace::propagator( const std::vector<double> &rmesh, double E, 
                         int l, double j ) {

    cmatrix_t ham = c_hamiltonian( rmesh, E, l, j ); 

    // Form of Propagator is 1 / ( E - H ), where H is the hamiltonian
    ham *= -1.0;
    cmatrix_t G( ham ); // Propagator

    for( unsigned int i = 0; i < rmesh.size(); ++i ) G( i, i ) += E;
        
    numeric::linalg::inverse_mtx( G );
    return G;
}

// This function is used in the two functions that occur after this one
bool 
compare_N( std::pair< double, std::vector<double> > pair1, 
           std::pair< double, std::vector<double> > pair2 ) {

    if( pair1.first < pair2.first ) return true;
    else return false;

}

// Find the real eigenvalues by diagonalizing Hamiltonian
// Calculates the eigenvalues only
std::vector< double >
boundRspace::real_eigvals( const matrix_t &Ham ) {

    // Find eigenvalues
    cvector_t evals = numeric::linalg::eigenvalues( Ham );
    std::list< double > real_evals;

    // Find the negative solutions
    for ( unsigned int m = 0; m < evals.size(); ++m ) {

        real_evals.push_back( real( evals( m ) ) );
    }
    
    // Check to make sure that eigenvalues were found
    if ( real_evals.empty() ) {
        std::cout << "No Real Eigenvalues Found!" << std::endl;
    }

    real_evals.sort();

    // return a vector instead of a list
    std::list< double >::iterator it;
    std::vector< double > evals_vec;
    for( it = real_evals.begin(); it != real_evals.end(); ++it ) {

        evals_vec.push_back( *it );
    }

    return evals_vec;

}

// Find the real states by diagonalizing Hamiltonian
// Calculates the eigenvalues and the eigenvectors
std::vector< eigen_t >
boundRspace::real_eigvecs( const matrix_t &Ham ) {

    // Find eigenvalues and eigenvectors
    std::pair< cvector_t, matrix_t > eval_evec = numeric::linalg::eig( Ham );
    cvector_t evals = eval_evec.first;
    matrix_t evecs = eval_evec.second;

    // find the eigenvalue corresponding to N
    std::list< eigen_t > real_list;

    for ( unsigned int m = 0; m < evals.size(); ++m ) {

        std::vector< double > eigvec;
        for ( unsigned int i = 0; i < evecs.size1(); ++i ) {

            eigvec.push_back( evecs( i, m ) );
        }
            
        real_list.push_back( std::make_pair( real( evals( m ) ), eigvec ) );
    }

    real_list.sort( compare_N );

    // return a vector instead of a list
    std::list< eigen_t >::iterator it;
    std::vector< eigen_t > real_vec;
    for ( it = real_list.begin(); it != real_list.end(); ++it ) {

        real_vec.push_back( *it );
        
    }

    return real_vec;

}

// Calculates the eigenvalues of the Schroedinger-like equation
// self-consistently
double
boundRspace::find_level( const std::vector<double> &rmesh, double Estart, 
                         int N, int l, double j, double tol ) {

    // construct Hamiltonian with starting energy value
    matrix_t ham = re_hamiltonian( rmesh, Estart, l, j );

    // get the negative eigenvalues sorted from lowest to highest
    std::vector< double > evals = real_eigvals( ham );

    std::vector<double> Ein_vec;
    double Ein = evals.at( N );

    ham.clear();
    evals.clear();
    ham = re_hamiltonian( rmesh, Ein, l, j );
    evals = real_eigvals( ham );

    double Eout = evals.at( N );

    // find self consistent eigenvalues for each principal quantum number N
    double tol2 = 2 * tol; 

    int iter = 0;
    int iter_max = 30;
    int iter_med = 8;

/*    std::cout << "------------------------------------------" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "n = " << N << " l = " << l << " j = " << j << std::endl;
    std::cout << " " << std::endl;
*/
    while ( std::abs( Eout - Ein ) > tol ) {

        iter += 1;
//        std::cout << iter << " " << Ein << " " << Eout << std::endl;

        // Assign new input energy
        if ( ( iter > iter_max ) && ( std::abs( Eout - Ein ) > tol2 ) ) {
            
            std::cout << "Eigenvalue not converging." << std::endl;
            break;
        }
        else if ( ( iter > iter_med ) && 
                  ( std::abs( Eout - Ein ) > tol2 ) ) {

            Ein = ( Eout + Ein ) / 2.;
        }
        else Ein = Eout;

        // construct Hamiltonian with new energy value
        ham.clear();
        ham = re_hamiltonian( rmesh, Ein, l, j );

        // get new eigenvalues
        std::vector< double > evals2 = real_eigvals( ham );
        if( evals2.empty() ) continue;

        Eout = evals2.at( N ); // New output energy

    } // end self-consistency loop

    ham.clear();
    ham = re_hamiltonian( rmesh, Eout, l, j );

    std::vector< double > evals3 = real_eigvals( ham );
    double eval = evals3.at( N );

    // check stability
    if ( std::abs( eval - Eout ) > tol ) {
  /*      std::cout << "Eigenvalue not converged within " << tol 
                  << " MeV." << std::endl;
        std::cout << "Difference is " << eval - Eout << " MeV." << std::endl;
   */ }

    return Eout;
}

// This function returns an overlap function this is normalized to 1. 
// The input eigenvector is actually the overlap function multiplied by r,
// so the eigenvector is divided by r in this function.
std::vector< double > 
boundRspace::normalize( const std::vector< double > &rmesh, 
           const std::vector<double> &eigenvector ) {

    double delt = rmesh.at(1) - rmesh.at(0);

    double norm2 = 0;
    for( unsigned int i = 0; i < rmesh.size(); ++i ) {

        norm2 += std::pow( eigenvector[ i ], 2 ) * delt;
        
    }
    double norm = std::sqrt( norm2 );
    std::vector< double > efx;
    double sign = 1;
    if ( eigenvector[0] < 0 ) sign = -1;
    for( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        efx.push_back( sign * eigenvector[i] / norm / rmesh[i] );
        //cout<<rmesh[i]<<" "<<efx[i]<<endl;
    }

    return efx;
}

eigen_t
boundRspace::find_boundstate( const std::vector<double> &rmesh, double Estart, 
                              int N, int l, double j, double tol ) {
     //   for (int kk=0;kk<rmesh.size();++kk){  
       // std::cout << rmesh[kk]<<std::endl;}
//	std::cout <<Estart<<" " <<N<<" "<<l<<" "<<j<<" " <<tol<<std::endl; 
        // get level
        double QPE = find_level( rmesh, Estart, N, l, j, tol );
      //  matrix_t ham = re_hamiltonian( rmesh, QPE, l, j );
        matrix_t ham = re_hamiltonian( rmesh, QPE, l, j );

        std::vector< eigen_t > eig_info = real_eigvecs( ham );

        double eval = eig_info[N].first;
        std::vector<double> &evec = eig_info[N].second;

        if ( std::abs( eval - QPE ) > tol ) {
            std::cout << "Functions 'numeric::linalg::eigenvalues' " 
                      << "and " <<  "'numeric::linalg::eig' " 
                      << "not " << "consistent with each other." << std::endl;
            
            std::cout << "From eigenvalues: " << QPE << std::endl;
            std::cout << "From eig: " << eval << std::endl;


            // choose level that is outside continuum
            ham.clear();
            ham = re_hamiltonian( rmesh, eval, l, j );
            std::vector< eigen_t > eig_info2 = real_eigvecs( ham );
            double eval2 = eig_info2[N].first;
            std::vector< double > &evec2 = eig_info2[N].second;

            if ( eval2 > eval ) {

                for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                    evec[i] = evec2[i];
                }
                eval = eval2;
            }

            std::cout << "New Calculation: " << eval << std::endl;
            //std::cout << "Emax: " << Emax << std::endl;

        }
        // Normalize Eigenfunction
        std::vector<double> efx = normalize( rmesh, evec );

        return std::make_pair( eval, efx );
}


// Calculate the Spectroscopic Factor. 'QPE' is the quasiparticle
// ( or quasihole ) energy and 'QPF' is the corresponding overlap function
double 
boundRspace::sfactor( const std::vector<double> &rmesh, double QPE, int l, 
                      double xj, const std::vector<double> &QPF ) {

    Pot->setEnergy( QPE );
    Pot->setAM( l, xj );

    double delt = rmesh.at(1) - rmesh.at(0);

    double dsigma_local = 0;
    double dsigma_nonlocal = 0;

    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        // local parts
        double dV_local = Pot->der_disp_localPart( rmesh[i] );

        dsigma_local += delt * std::pow( rmesh[i], 2 ) * dV_local 
                      * std::pow( QPF[i], 2 );
                      
        // nonlocal parts
        for ( unsigned int k = 0; k < rmesh.size(); ++k ) {
            
            double dV_nonlocal = 
                Pot->der_disp_nonlocalPart( rmesh[i], rmesh[k] );

            // nonlocal part is already multiplied by rr'
            dsigma_nonlocal += delt * delt * dV_nonlocal
                             * rmesh[i] * QPF[i] * rmesh[k] * QPF[k];
        }
    }

    double dsigma = dsigma_local + dsigma_nonlocal;
    return 1.0 / ( 1.0 - dsigma );
}

double 
boundRspace::rms_radius( const std::vector<double> &rmesh, 
                         const std::vector<double> &QPF ) {

    double delt = rmesh.at(1) - rmesh.at(0);

    double rsq = 0;

    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        rsq += std::pow( rmesh[i], 4 ) * delt * std::pow( QPF[i], 2 );
    }

    return std::sqrt( rsq );
}

// Occupation calculated from the density matrix. 'd_mtx' is the 
// density matrix for a specific value of L and J
double 
boundRspace::occupation( const std::vector<double> &rmesh, 
                         const matrix_t &d_mtx, 
                         const std::vector<double> &QPF ) {


    double delt = rmesh.at(1) - rmesh.at(0);

    double occ = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        for( unsigned int j = 0; j < rmesh.size(); ++j ) {

            occ += std::pow( rmesh[i], 2 ) * QPF[i] * d_mtx( i, j )
                 * std::pow( rmesh[j], 2 ) * QPF[j] * std::pow( delt, 2 );
        }
    }

    return occ;
}

double
boundRspace::spectral_function_k_space( const std::vector< double > &rmesh, 
                                        double E, double k ) {

    double delt = rmesh.at( 1 ) - rmesh.at( 0 );

    double S_of_kE = 0;
    for( int L = 0; L < Lmax + 1; ++L ) {

        for ( int up = -1; up < 2; up += 2 ) {

            double J = L + up / 2.0;
            if ( J < 0 ) continue;

            // Propagator
            cmatrix_t G = propagator( rmesh, E, L, J );

            // Spectral Function in momentum space, S( k; E ) 
            double rsums = 0;
            for( unsigned int i = 0; i < rmesh.size(); ++i ) {

                double rho1 = k * rmesh[i];
                double jl1 = gsl_sf_bessel_jl( L, rho1 );

                for( unsigned int j = 0; j < rmesh.size(); ++j ) {

                    double rho2 = k * rmesh[j];
                    double jl2 = gsl_sf_bessel_jl( L, rho2 );
                    
                    rsums -= rmesh[i] * jl1 * imag( G( i, j ) ) 
                           * rmesh[j] * jl2 * delt * 2 / M_PI / M_PI;
                }

            } // end loop over radial coordinates

            S_of_kE += ( 2 * J + 1 ) * rsums;
            
        } // end loop over J
    } // end loop over L

    return S_of_kE;
}

// a simple routine to find the absolute maximum of a function
std::pair< double, double >
boundRspace::find_max( const std::vector< double > &x_values, 
                       const std::vector< double > &y_values ) {

    if ( x_values.size() != y_values.size() ) {

        std::cout << "In function find_max: " 
                  << "x_values and y_values have different sizes!" 
                  << std::endl;

        std::cout << "x_values.size() = " << x_values.size() << std::endl;
        std::cout << "y_values.size() = " << y_values.size() << std::endl;

        std::abort();
    }

    double y_coordinate = y_values.front();

    unsigned int i_at_maximum = 0;
    for ( unsigned int i = 0; i < x_values.size(); ++i ) {

        if ( y_values[i] > y_coordinate ) {

            y_coordinate = y_values[i];
            i_at_maximum = i;
        }
    }

    return std::make_pair( x_values[i_at_maximum], y_values[i_at_maximum] );
}

// rough estimate of the width of a quasi-hole state
double 
boundRspace::approx_width( const std::vector<double> &rmesh, double QPE, 
                           int L, double J, const std::vector< double > &QPF ) {

    Pot->setEnergy( QPE );
    Pot->setAM( L, J );

    double delt = rmesh.at( 1 ) - rmesh.at( 0 );

    double sum = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        for ( unsigned int j = 0; j < rmesh.size(); ++j ) {

            // potential is already multiplied by rr'
            sum += delt * delt * rmesh[i] * rmesh[j] 
                 * QPF[i] * QPF[j] * Pot->nonlocalIM( rmesh[i], rmesh[j] );
        }
    }

    return 2 * std::abs( sum );
}

// a simple routine to calculate the width of quasi-hole state
// from the spectral strength function, s_of_E_values, which should
// be folded with the relevant wave function so that there is only
// one peak. 
double
boundRspace::find_width( const std::vector< double > &e_values,
                         const std::vector< double > &s_of_E_values ) {

    std::pair< double, double > qh_peak = find_max( e_values, s_of_E_values );
    double half_max_y = qh_peak.second / 2;

    double half_max_x = qh_peak.first - 10;
    for ( unsigned int i = 0; i < s_of_E_values.size(); ++i ) {

        if ( s_of_E_values[i] > half_max_y ) {

           half_max_x = ( e_values[i] + e_values.at( i - 1 ) ) / 2.0; 

           break;
        }
    }

    return 2 * std::abs( qh_peak.first - half_max_x );
    
}

// Function for storing the quasi-hole levels
// lj_eigen_t type defined in boundRspace.h
std::vector< lj_eigen_t > 
boundRspace::get_bound_levels( const std::vector< double > &rmesh, 
                               double tol ) {

    std::cout << "Finding Levels" << std::endl;

    std::vector< lj_eigen_t > lj_levels;
    // loop over L
    for ( int L = 0; L < Lmax + 1; ++L ) {

        // loop over J
        for ( int up = -1; up < 2; up +=2 ) {

            double J = L + up / 2.0;
            if ( J < 0 ) continue;

            std::vector< eigen_t > levels;
            for ( int N = 0; N < max_n; ++N ) {

                double Estart = Efermi; // for now
                eigen_t bound_info = 
                    find_boundstate( rmesh, Estart, N, L, J, tol );

            //    cout << N << " " << L << " " << "  " << J << " " 
            //       << bound_info.first << endl ;     

                if ( bound_info.first >= 0 ) break; 
                levels.push_back( bound_info );
                    
            } // End loop over N

            lj_levels.push_back( levels );

        } // End loop over J
    } // End loop over L

    return lj_levels;

}

// function for creating the energy meshes for each LJ. The meshes depend
// on whether there are bound states and how close a bound state is to 
// the energy at which the imaginary part begins.
//
// 'rmesh' contains the radial grid
// 'tol' specifies the accuracy to which the bound states should be calculated
//
// type 'mesh_t' is defined in the boundRspace header file

std::vector< mesh_t >
boundRspace::get_emeshes( const std::vector< double > &rmesh, 
                          double Emin, double Emax, 
                          const std::vector< lj_eigen_t > &bound_levels ) {

    std::cout << "Creating energy meshes ... " << std::endl;

    // vector for storing the energy meshes for each LJ combination
    std::vector< mesh_t > emesh_vec;

    // loop over L
    for ( int L = 0; L < Lmax + 1; ++L ) {

        // loop over J
        for ( int up = -1; up < 2; up +=2 ) {

            double J = L + up / 2.0;
            if ( J < 0 ) continue;

            int lj_index = index_from_LJ( L, J );
            const std::vector< eigen_t > &levels = bound_levels.at( lj_index );

            lj_eigen_t qh_in_continuum_vec;
            for ( unsigned int N = 0; N < levels.size(); ++N ) {

                if ( levels[N].first >= 0 ) { 
                    
                    std::cout << "Error: level >= 0." << std::endl;
                    std::cout << "N, L, J = " << N << " " << L << " " << J
                              << std::endl;

                    std::cout << "Level = " << levels[N].first << std::endl;
                    
                    std::abort();
                }

                // the locations of the quasihole peaks that
                // are in or near the continuum will be used
                // to construct the energy mesh. 
                if ( levels[N].first < ( Emax + 3 ) ) {
                    
                    qh_in_continuum_vec.push_back( levels[N] );

                    if ( levels[N].first < Emax ) {
//
                    std::cout << "N, L, J = " << N << " " << L << " " << J
                              << std::endl;
                        std::cout << " " << std::endl;
                        std::cout << "Level in continuum." << std::endl;
                        std::cout << " " << std::endl;
                    }
                    else {

                        std::cout << " " << std::endl;
                        std::cout << "Level outside but near continuum." 
                                 << std::endl;
                       std::cout << " " << std::endl;
                    }
                }

            } // End loop over N

            // Create Energy Grid for calculating density matrix etc.
            // based on the location of the quasihole peaks
            double Elowmax = -60; // default
            double Emedmax = Emax - ph_gap; // default
            double Edeltlow =5; // default
            double Edeltmed =1;// default
            double Edelthigh =.2; // default

            if ( qh_in_continuum_vec.size() == 1 ) {
                double QPE = qh_in_continuum_vec[0].first;
                std::vector< double > &QPF = qh_in_continuum_vec[0].second;

                // get rough estimate of quasihole width
                double gamma = approx_width( rmesh, QPE, L, J, QPF );

                if ( QPE > Emedmax ) {
                    
                    Edelthigh =0.05;
                //    Edelthigh =0.001;
                }
                else {

                    Elowmax = QPE - gamma;
                 //  Edeltmed =  0.001;
                   Edeltmed =  0.2;
                }

            }
            else if ( qh_in_continuum_vec.size() > 1 ) {
            cout<<"HEY"<<endl;
                // Find the minimum and maximum bound-state
                // energies located in or near the continuum
                eigen_t &bound_max = qh_in_continuum_vec.back();
                eigen_t &bound_min = qh_in_continuum_vec.front();

                double QPE_max = bound_max.first;
                double QPE_min = bound_min.first;

                if ( QPE_max > Emedmax ) 
			Edelthigh =0.05;
	//		Edelthigh =0.001;

                // get rough estimate of quasihole width of 
                // lowest level
                double gamma = 
                    approx_width( rmesh, QPE_min, L, J, bound_min.second );

                Elowmax = QPE_min - gamma;
               Edeltmed = 0.2 ;
            //    Edeltmed = 0.02 ;

            }

            mesh_t emesh = energy_mesh_fast( Emin, Elowmax, Emedmax, Emax, 
                                             Edeltlow, Edeltmed, Edelthigh );

            emesh_vec.push_back( emesh );

        } // End loop over J
    } // End loop over L

    if ( emesh_vec.size() != static_cast< unsigned >( 2 * Lmax + 1 ) ) {

        std::cout << "Some  LJ combinations are missing." << std::endl;
        std::cout << "emesh_vec.size() = " << emesh_vec.size() << std::endl;
        std::cout << "expected size = " << 2 * Lmax + 1 << std::endl;
    }

    std::cout << "Finished creating energy meshes." << std::endl;

    return emesh_vec;
}

int
boundRspace::L_from_index( unsigned int index ) {

    int L;
    if ( index % 2 == 0 ) L = static_cast<int>( index / 2 );
    else L = static_cast<int>( ( index + 1 ) / 2 );

    return L;
}

double
boundRspace::J_from_index( unsigned int index ) {

    double J;
    if ( index % 2 == 0 ) J = static_cast< double >( index + 1 ) * 0.5;
    else J = static_cast< double >( index ) * 0.5;

    return J;
}

int
boundRspace::index_from_LJ( int L, double J ) {

    // check to make sure that L is not greater than Lmax
    if ( L > Lmax ) {

        std::cout << "in function 'index_from_LJ': "
                  << "L too large." << std::endl;

        std::cout << "L = " << L << std::endl;
        std::cout << "Lmax = " << Lmax << std::endl;
    }

    int index;
    if ( J > static_cast<double>( L ) ) index = 2 * L;
    else index = 2 * L - 1;

    return index;

}

mesh_t
boundRspace::get_lj_emesh( int L, double J, 
                           const std::vector< mesh_t > &emesh_vec ) {

    int lj_index = index_from_LJ( L, J );
    return emesh_vec.at( lj_index );
}

std::vector< double > 
boundRspace::get_lj_energy_vector( int L, double J, 
                                   const std::vector< mesh_t > &emesh_vec ) {

    int lj_index = index_from_LJ( L, J );

    std::vector< double > vec;
    for ( unsigned int n = 0; n < emesh_vec.at( lj_index ).size(); ++n ) {
        
        vec.push_back( emesh_vec[lj_index][n].first );

    }

    return vec;
}


std::vector< prop_t > // type defined in header file
boundRspace::get_propagators( const std::vector< double > &rmesh, 
                              const std::vector< mesh_t > &emesh_vec ) {

    std::cout << "Calculating propagators ..." << std::endl;

    std::vector< prop_t > prop_vec;
    for ( unsigned int m = 0; m < emesh_vec.size(); ++m ) {
        
        const mesh_t &emesh = emesh_vec[m];

        int L = L_from_index( m );
        double J = J_from_index( m );

        prop_t G_vec;
        for ( unsigned int n = 0; n < emesh.size(); ++n ) {

            double E = emesh[n].first;

//              cout<<E <<endl;
            // Store propagator for different energy values in vector
            G_vec.push_back( propagator( rmesh, E, L, J ) );
/*
// ***************NEW to test*************************************************************
            if(n==10 && m==0){ cmatrix_t  prop_test = propagator (rmesh,E,L,J);
				for (int ii=0 ; ii< rmesh.size(); ii++) { 
					for (int jj=0; jj<rmesh.size(); jj++) { 
						if (ii==1 && jj==1) {std::cout <<E <<" "<< m << " " << prop_test(ii,jj) << std::endl;}
					 }
				}
	    }	
// ****************************************************************************										
*/
        } // end loop over energy

   // std::cout << "Calculating propagators2 , end loop over energy ..." << std::endl;
        // Store propagator for different lj in vector
        prop_vec.push_back( G_vec );
    }

    std::cout << "Finished calculating propagators." << std::endl;
    std::cout << " " << std::endl;

    return prop_vec;
}

// calculates the spectral strength S( E ) from the propagator
std::vector< double >
boundRspace::spectral_strength( int L, double J,
                                const std::vector< double > &rmesh, 
                                const std::vector< mesh_t > &emesh_vec, 
                                const std::vector< prop_t > &prop_vec ) {

    // select energy mesh and propagator specified by L and J
    int lj_index = index_from_LJ( L, J );
    const mesh_t &emesh = emesh_vec.at( lj_index );
    const prop_t &propE = prop_vec.at( lj_index );


    std::vector< double > spf;
    for ( unsigned int n = 0; n < emesh.size(); ++n ) {

        // Spectral Function 
        complex_t trace = 0;
        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                    
            trace += propE[n]( i, i ); 
        }

        spf.push_back( -imag( trace ) / M_PI ); 
        
    }

    return spf;

}

// calculates the spectral strength dS(r, E )/dr from the propagator

std::vector< std::vector< double > >
boundRspace::spectral_strength_der( int L, double J,
                                const std::vector< double > &rmesh, 
                                const std::vector< mesh_t > &emesh_vec, 
                                const std::vector< prop_t > &prop_vec ) {

    // select energy mesh and propagator specified by L and J
    int lj_index = index_from_LJ( L, J );
    const mesh_t &emesh = emesh_vec.at( lj_index );
    const prop_t &propE = prop_vec.at( lj_index );
    std::vector< std::vector< double> > spec_derivative; 
    std::vector< double > trac;
    for ( unsigned int n = 0; n < emesh.size(); ++n ) {

        // derivative of Spectral Function dS(r,e)/dr 
        for( unsigned int i = 0; i < rmesh.size()-2; ++i ) {
          
	   trac.push_back(imag((propE[n](i+1,i+1)/((rmesh[i+2]-rmesh[i+1]) * rmesh[i+1] * rmesh[i+1])- propE[n](i,i)/(rmesh[i] * rmesh[i] * (rmesh[i+1]-rmesh[i])))/(rmesh[i+1]-rmesh[i]))); 	 
          }
	    spec_derivative.push_back(trac);
    }
     
   return spec_derivative;

}

// calculates the spectral strength S( E ) from the propagator
// folded with a quasi-hole wavefunction
std::vector< double >
boundRspace::spectral_strength_QH( int L, double J, 
                                   const std::vector< double > &rmesh, 
                                   const std::vector< mesh_t > &emesh_vec, 
                                   const std::vector< prop_t > &prop_vec,
                                   const std::vector< double > &QHF ) {

    // select energy mesh and propagator specified by L and J
    int lj_index = index_from_LJ( L, J );

    const mesh_t &emesh = emesh_vec.at( lj_index );
    const prop_t &propE = prop_vec.at( lj_index );

    std::vector< double > spf_qh;
    for ( unsigned int n = 0; n < emesh.size(); ++n ) {

        double fspf = 0;
        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
        for( unsigned int j = 0; j < rmesh.size(); ++j ) {

            fspf -= rmesh[i] * rmesh[j] * rdelt / M_PI
                  * QHF[i] * QHF[j] * imag( propE[n]( i, j ) );

        } // end loop over i
        } // end loop over j

        spf_qh.push_back( fspf ); 
        
    }

    return spf_qh;

}

double 
boundRspace::particle_number( const std::vector< double > &rmesh, double Emax,
                              const std::vector< mesh_t > &emesh_vec,
                              const std::vector< prop_t > &prop_vec,
                              const std::vector< lj_eigen_t > &bound_levels ) {

    std::ofstream file( "Output_test/strengths.out" );
    // loop over L and J
    double number = 0;
    for ( unsigned int i = 0; i < emesh_vec.size(); ++i ) {

        const mesh_t &emesh = emesh_vec.at( i );
        const prop_t &G = prop_vec.at( i );
        const lj_eigen_t &levels = bound_levels.at( i );

        int L = L_from_index( i );
        double J = J_from_index( i );

        // contribution from continuum 
        double lj_strength = 0;
        for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                
            double Edelt = emesh[m].second;

            // Spectral Function 
            complex_t trace = 0;
            for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                    
                trace += G[m]( i, i ); 
            }

            double spf = -imag( trace ) / M_PI; 
            lj_strength += Edelt * spf;

        } // end loop over energy
            
        //  contribution from quasi-hole peaks outside of the continuum 
        for ( unsigned int N = 0; N < levels.size(); ++N ) {

            if ( ( levels[N].first <= Efermi ) && 
                 ( levels[N].first > Emax ) ) {
                
                // Quasiparticle energy
                double QPE = levels[N].first; 

                // Quasiparticle wavefunction
                const std::vector<double> &QPF = levels[N].second;

                // Spectroscopic Factor
                double S = sfactor( rmesh, QPE, L, J, QPF );

                lj_strength += S;

            } // endif

        } // End loop over N

        number += ( 2 * J + 1 ) * lj_strength;

        file << "Strength for LJ = " << L << " " << J 
             << ": " << lj_strength << " " << number << std::endl;

    } // end loop over L and J

    file.close();
    file.clear();

    return number;

}

std::vector< double >
boundRspace::charge_density( const std::vector< double > &rmesh, double Emax,
                             const std::vector< mesh_t > &emesh_vec,
                             const std::vector< prop_t > &prop_vec,
                             const std::vector< lj_eigen_t > &bound_levels ) {

    std::vector<double> point_dist;
    point_dist.assign( rmesh.size(), 0 );

    // loop over L and J
    for ( unsigned int i = 0; i < emesh_vec.size(); ++i ) {

        const mesh_t &emesh = emesh_vec.at( i );
        const prop_t &G = prop_vec.at( i );
        const lj_eigen_t &levels = bound_levels.at( i );

        int L = L_from_index( i );
        double J = J_from_index( i ); 
        // Loop over energy
        matrix_t d_mtx( rmesh.size(), rmesh.size() ); // density matrix
        d_mtx.clear();
        for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                
            double Edelt = emesh[m].second;

            // Density Matrix
            for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                    
                    d_mtx( i, j ) -= Edelt * imag( G[m]( i, j ) ) / M_PI
                                   / ( rmesh[i] * rmesh[j] * rdelt );
                }
            }


        } // end loop over energy
            
        // Charge Density from continuum
        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
    
            point_dist[i] += ( 2 * J + 1 ) * d_mtx( i, i ) / ( 4 * M_PI );
        }

        // Charge density from quasi-hole peaks outside of the continuum 
        for ( unsigned int N = 0; N < levels.size(); ++N ) {

            if ( ( levels[N].first <= Efermi ) && 
                 ( levels[N].first > Emax ) ) {
                
                // Quasiparticle energy
                double QPE = levels[N].first; 

                // Quasiparticle wavefunction
                const std::vector<double> &QPF = levels[N].second;

                // Spectroscopic Factor
                double S = sfactor( rmesh, QPE, L, J, QPF );

                for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                    point_dist[i] += ( 2 * J + 1 ) * S * QPF[i] * QPF[i] 
                                   / ( 4 * M_PI );

                } // end loop over r

                std::cout << "added to charge density " << N << " " 
                          << L << " " << J << " " << QPE << " " 
                          << S << std::endl;

            } // endif

        } // End loop over N

    } // end loop over L and J

    // Find normalization of point distribution
    double point_norm = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        point_norm += 4 * M_PI * point_dist[i] * std::pow( rmesh[i], 2 )
                    * rdelt;
    }

    double norm_fac = Z / point_norm;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    
        point_dist[i] *= norm_fac;
    }

    std::cout << "Folding density distribution" << std::endl;

    // Fold point charge density with neutron and proton charge distribution
    // FIXME can only use the same point distribution for neutrons and protons
    // if the nucleus has N = Z. Eventually, the code below needs to be made
    // more general.

    std::vector< double > rweights;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        rweights.push_back( rdelt );
    }

    std::vector<double> proton_folding = 
        folded_ch_density( rmesh, rweights, point_dist, 0.5, A );

    std::vector<double> neutron_folding = 
        folded_ch_density( rmesh, rweights, point_dist, -0.5, A );

    std::vector< double > chdf; // folded charge density
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        chdf.push_back( proton_folding[i] + neutron_folding[i] );
    }

    // FIXME still need to take into account spin-orbit correction
    return chdf;

}

// Calculates the point density distribution without normalizing to
// the experimental particle number
std::vector< double >
boundRspace::point_distribution(const std::vector< double > &rmesh, 
                                 double Emax,
                                 const std::vector< mesh_t > &emesh_vec,
                                 const std::vector< prop_t > &prop_vec,
                                const std::vector< lj_eigen_t > &bound_levels) {

    std::vector<double> point_dist;
    point_dist.assign( rmesh.size(), 0 );

     cout << "Emax = " << Emax <<"  " << emesh_vec.size() << endl ;

    // loop over L and J
    for ( unsigned int i = 0; i < emesh_vec.size(); ++i ) {

        const mesh_t &emesh = emesh_vec.at( i );
        const prop_t &G = prop_vec.at( i );
        const lj_eigen_t &levels = bound_levels.at( i );

        int L = L_from_index( i );
        double J = J_from_index( i );
        // Loop over energy
        matrix_t d_mtx( rmesh.size(), rmesh.size() ); // density matrix
        d_mtx.clear();
        for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                
            double Edelt = emesh[m].second;


            // Density Matrix
            for( unsigned int ii = 0; ii < rmesh.size(); ++ii ) {
                for( unsigned int j = 0; j < rmesh.size(); ++j ) {
//		   if ( (L >= 3  && emesh[m].first >= -14.6)){
//		 	d_mtx( ii,j) +=0;
//	    }
//		   else  { 
                    d_mtx( ii, j ) -= Edelt * imag( G[m]( ii, j ) ) / M_PI
                               / ( rmesh[ii] * rmesh[j] * rdelt );
//		   }
		}
            }


        } // end loop over energy
        // Charge Density from continuum
        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
    
            point_dist[i] += ( 2 * J + 1 ) * d_mtx( i, i ) / ( 4 * M_PI );

        }

        // Charge density from quasi-hole peaks outside of the continuum 
        for ( unsigned int N = 0; N < levels.size(); ++N ) {
      

            if ( ( levels[N].first <= Efermi ) && 
                 ( levels[N].first > Emax ) ) {
                
                // Quasiparticle energy
                double QPE = levels[N].first; 

                // Quasiparticle wavefunction
                const std::vector<double> &QPF = levels[N].second;

                // Spectroscopic Factor
                double S = sfactor( rmesh, QPE, L, J, QPF );
                for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                    point_dist[i] += ( 2 * J + 1 ) * S * QPF[i] * QPF[i] 
                                   / ( 4 * M_PI );

                } // end loop over r
                
                std::cout << "added to charge density " << N << " " 
                          << L << " " << J << " " << QPE << " " 
                          << S <<" " <<std::endl;
            } // endif

        } // End loop over N

    } // end loop over L and J

    return point_dist;

}

double
boundRspace::chd_rms_radius( const std::vector< double > &rmesh,
                             const std::vector< double > &chd ) {

    // Calculate RMS Radius
 //   double norm = 0; // norm 
    double norm = 0.; // norm 
    double sum_r2 = 0; // integral of charge distribution weighted with r * r
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        norm += std::pow( rmesh[i], 2 ) * chd[i] * rdelt;
        sum_r2 += std::pow( rmesh[i], 4 ) * chd[i] * rdelt;
//        if (rmesh[i] > 7.0) break;

    }

    return std::sqrt( sum_r2 / norm );

}

