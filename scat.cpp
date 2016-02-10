#include "scat.h"
#include "waves.h"




//***********************************************************************
double const scat::kconstant =  .048192;
double const scat::e2 = 1.44; // Coulomb constant
double const scat::m0 = 931.5; // nucleon mass
double const scat::pi = acos(-1.);
double const scat::EminRel = 1.;
int const scat::Nnumerov = 200;
int const scat::Imethod = 1;
bool const scat::local = false;
//*************************************************************************
//default constructor
scat::scat(){}
//**************************************************************************
  /**
   * loads the potential so subsequent operations
\param HartreeFock points to the Hartree-Fock potential
\param Volume0 point to the volume imaginary potential
\param Surface0 points to the  Surface imaginary pot
\param SpinOrbit0 points to the real and imaginary spin-orbit potential
\param Rc0 is the radius in fm for the Coulomb potential
\param energCM0 is the energy in the center-of-mass frame in units of MeV
\param energyLab0 is the energy in the labortory frame
   */
void scat::loadPotential(hartreeFock *HartreeFock0, 
	   volume *Volume0, surfaceTF *Surface0, 
           spinOrbit * SpinOrbit0, double Rc0,double energyCM0, 
			 double energyLab0,double beta_nl_R00, 
                         double beta_nl_R10,double beta_nl_I0)

{
  beta_nl_R0 = beta_nl_R00;
  beta_nl_R1 = beta_nl_R10;
  beta_nl_I = beta_nl_I0;
  HartreeFock = HartreeFock0;
  if (HartreeFock->A != 0.)beta_max = max(beta_nl_R0,beta_nl_R1);
  else beta_max = beta_nl_R0;
  beta_max = max(beta_max,beta_nl_I);

  Volume = Volume0;
  Surface = Surface0;
  SpinOrbit = SpinOrbit0;
  Rc = Rc0;

  energyCM = energyCM0;
  energyLab = energyLab0;
  Kwave2 = kconstant*mu*energyCM;


  CoulDis = 1.73*Z*Zp/Rc; // Coulomb Displacement

  //relativistic correction
  if (energyCM > EminRel) Kwave2 = kconstant*pow(A/(1.+ A + energyCM/m0),2)*
		       energyLab*(energyLab/2./m0 + 1.);  

  Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1

  if (energyCM > EminRel) gammaRel = 2.*(energyCM/m0 +1.)/(energyCM/m0 + 2.);
  else gammaRel = 1.;

  if (energyCM > 0.) muhbar = Kwave2/energyCM;
  else muhbar = kconstant*mu;

  gamma = fabs(gammaRel*Z*Zp*e2*Kwave/energyCM/2.); //Sommerfeld parameter
  //gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter
  konst = pi/pow(Kwave,2); //const for cross sections in units of Fermi squared
  konst *= 10.; //umits of mb

}
//***************************************************************************
  /**
   * Constructor
\param Zp0 is atomic number of projectile(0=neutron,1=proton)
\param Z0 is atomic number of target
\param A0 is the mass number of the target
\param flag0 indicates integrating of wavefunctions takes place
  */
scat::scat(double Zp0, double Z0, double A0, bool flag0, string *title0)
  :compound(title0)
{
  init(Zp0,Z0,A0,flag0);
}
//************************************************************************
  /**
   * initializes for a new calculation
\param Zp0 is atomic number of projectile(0=neutron,1=proton)
\param Z0 is atomic number of target
\param A0 is the mass number of the target
\param flag0 indicates integrating of wave functions takes place
  */
void  scat::init(double Zp0,double Z0, double A0, bool flag0)

{
  Zp = Zp0;
  flag = flag0;
  A = A0;
  Z = Z0;
  mu = A/(1.+A);
  if (Zp!=0.)proton = 1;
  else proton = 0;
  if (flag == 0) return;
  elastic = 1;
  initIntegration();
}
//*************************************************************************
  /**
   *returns the real total potential (centripetal not included) in MeV
   *Includes the dispersive correction
   *It is scaled by gamma \f$\sqrt{1-\left(\frac{v}{c}\right)^2}\f$
\param r is the radial distance in fm
   */
double scat::RealPotential(double r)
{
  // first nuclear potential

  double one = HartreeFock->potential(r);
  double two = Volume->DispersiveCorrection(r);
  double three = Surface->DispersiveCorrection(r);
  double fact = one + two + three;



  if (proton == 1) // finally the Coulomb potential
    {
      if (r > Rc) fact += e2*Z*Zp/r;
      else fact += e2*Z*Zp/2./Rc*(3.-pow(r/Rc,2));
    }

  if (r > 0.) 
    {
      // next the spin-orbit potential
     fact += SpinOrbit->RealPotential(r);
    }

  return fact*gammaRel; 
}
//**************************************************************************
  /**
   * return the imaginary potental in MeV
   *It is scaled by gamma \f$\sqrt{1-\left(\frac{v}{c}\right)^2}\f$
   \param r is the radial distance in fm
   */
double scat::ImaginaryPotential(double r)
{
  double fact =  Surface->ImaginaryPot(r) 
     + Volume->ImaginaryPot(r);
  if (r > 0.) 
    {
      // next the spin-orbit potential
      fact += SpinOrbit->ImaginaryPotential(r);
    }
  return fact*gammaRel; 
}

//*****************************************************************
valarray<double> scat::diff(double r , valarray<double> w)
{
 //       REAL (kind=8) :: pot_real,pot_imag   
 // this subroutine is used by the mersion subroutine 
 // w[0] = REAL(u(r)) real part of  radial wave function
 //w[1] = REAL(du/dr) derrivative of real part of radial wave function
 //w[2] = IMAG(u(r)) Imaginary part of radial wave function
 //w[3] = IMAG(du/dr) derrivative of imaginary part
 //F(0:3) are the derivaties of w(0:3)




  int n = w.size();
  valarray<double> f(n);
  f[0] = w[1]; //these are equal by definition

  double potReal = RealPotential(r)*muhbar;
  if ( l > 0) potReal += (double)(l*(l+1))/pow(r,2);

  f[1] = -(Kwave2 - potReal)*w[0]; 


  if (n == 4)
    {
      f[2] = w[3]; // equal by definition
      double potImag = ImaginaryPotential(r)*muhbar;
      f[1] -= potImag*w[2]; 
      f[3] = -(Kwave2 - potReal)*w[2] + potImag*w[0];

    }

  return f;
}

//*****************************************************************
  /**
   * For positive energies, this Integrates the wave functions out from 
   * the origin to the matching radius, where the wavefunction is match to
   * Coulomb wavefunctions or spherical Bessel functions.
   * It determines the scattering phase shift, uses Merson Method 
   * to integrate wave function
   */

int scat::integrateWaveM()
{


  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave; 


  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];

  valarray<double> WaveFunct(4);

  l = 0;
  lStop = -1;
  for (;;) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;
          LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
          SpinOrbit->setAM(l,j);

          // potential at start
          double Vstart = RealPotential(rStart);
          double Wstart = ImaginaryPotential(rStart);
          //derivative of potential at start
          double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;
          double dWstart = (ImaginaryPotential(rStart+.01) - Wstart)/0.01;

	  // initialize wavefunctions
          double fact = pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[0] = pow(rStart,l+1) 
                                - muhbar*(energyCM-Vstart)*fact; // real part
          WaveFunct[2] =  Wstart*muhbar*fact;              //imaginary part

          // derivative of wavefunction
          fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
          WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                           - muhbar*(energyCM-Vstart)*fact; // real
          WaveFunct[3] = muhbar*Wstart*fact;          // imaginary

          fact = muhbar*pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[1] += dVstart*fact;
          WaveFunct[3] += dWstart*fact;


	  //integrate wavefunctions out to matching radius
          solveMerson(&WaveFunct,rStart,rStop);

          if (ok == 0) 
	    {
              cout << "j= " << j << " l= " << l << " Ecm= " <<
		     energyCM << endl;
              return 0;
	    }

          //outWave gives derivates with respect to rho = Kwave*r
	  //but at this point WaveFunctOut have derivative with respect to r
          // so make them with respect to rho=kwave*r

	  WaveFunct[1] /= Kwave;
          WaveFunct[3] /= Kwave;


	  // match wave functions 
	  //real WaveFunct = AA*F + BB*G
	  double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
	  double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

	  // imaginary part => Wavefunct  = CC*F + DD*G
	  double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
	  double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];


	  double denominator = pow(AA+DD,2) + pow(CC-BB,2);
	  double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double etaImag = 2.*(AA*BB+CC*DD)/denominator;

          int up = (plusMinus+1)/2;  //=0 form down and =1 for up
	  phaseShift[l][up] = atan2(etaImag,etaReal);
          if (phaseShift[l][up] < 0.) phaseShift[l][up] += 2.*pi;
          phaseShift[l][up] /=2.;
         
	  eta2[l][up] = pow(etaReal,2) + pow(etaImag,2);
          eta[l][up] = complex<double>(etaReal,etaImag);


	}
      //check to see if we have included enough l-waves
      if (1.-eta2[l][1] < 0.001 && 1.-eta2[l][0] < 0.001 )
	{
          lStop = l;
	  break;
	}
      if (isnan(eta2[l][1]) || isnan(eta2[l][0]))
	 {
	   lStop = l-1;
           //cout << 1. - eta2[lStop][1] << " " << 1. - eta2[lStop][0] << endl;
	   break;
	 } 
      l++;
      if (l > lMax) break;
    }
      if (lStop == -1) cout << "increase lMax" << endl;
      return 1;
}
//********************************************************
  /**
   * returns the elastic center-of-mass differntial cross section.
   * the function integrateWave() must be executed first
   \param theta center-of-mass angle in radians 
  */
double scat::DifferentialXsection(double theta)
{
  complex<double> A(0.,0.);
  complex<double> B(0.,0.);
  complex<double>tempA;
  complex<double>tempB;

  int ll = lMax;
  legendre Poly(ll);


  // start with l=0
  // eta is the S matrix

  A = eta[0][1];
  if (elastic) A -= 1.;

  if (proton) A *=  exp(complex<double>(0.,2.*Sigma[0])); 
  A *= Poly.LegendreP0(0,theta);


  //now the other l's

  for (int i=1;i<=lStop;i++)
    {
      //spin nonflip
      l = i;
      tempA = (double)(l+1)*eta[l][1] + (double)l*eta[l][0] ;
      if (elastic) tempA -=  complex<double>(2*l+1);
      tempA *= Poly.LegendreP0(l,theta);
      

      //spin flip
      tempB = (eta[l][1]-eta[l][0])*Poly.LegendreP1(l,theta);

      if (proton) // include Coulomb phase shift
	{
	  tempA *= exp(complex<double>(0.,2.*Sigma[l]));
          tempB *= exp(complex<double>(0.,2.*Sigma[l]));
	}

      A += tempA;
      B += tempB;
    }


  if (proton == 0 && theta < .01) cout << A << " " << B << endl;

  A /= complex<double>(0.,2.*Kwave);
  B /= complex<double>(0.,2.*Kwave);  



  if (proton)
    {
      complex<double>Acoulomb = -gamma/2./Kwave/pow(sin(theta/2.),2)*
	exp(complex<double>(0.,2.*Sigma[0] 
        - 2.*gamma*log(sin(theta/2.))));
      A += Acoulomb;
    }

  //analysing power
  complex<double>temp = 2.*A*conj(B);
  AnalyzePower = -imag(temp)/(pow(abs(A),2)+pow(abs(B),2));

  SpinRotation = real(temp)/(pow(abs(A),2)+pow(abs(B),2));



  double xsec = pow(abs(A),2) + pow(abs(B),2); //xsec in fm**2/sr
  //return xsec/100.; // units barns/sr
  return xsec*10.;    // units of mb/sr
}
//*****************************************************************
  /**
   * returns the absorbtion cross section in mb
   * also calculates the array scat.SigmaAb - i.e. the cross section as 
   * a function of compound nucleus spin
   * the function integrateWave() must be executed first
   */
double scat::AbsorptionXsection()
{
  //l==0 contribution
  SigmaAb[0] = (1. - eta2[0][1]) *konst;
  SigmaAb_pos[0] = SigmaAb[0];
  SigmaAb_neg[0] = 0.; 
  double tot = SigmaAb[0];
  //contribution from higher l waves
 
  //cout << 0 << " " << (1.-eta2[0][1]) << endl;

  for (int i=1;i<=lStop;i++)
    {
      l = i;
      int parity = 1-2*(i%2);
      double up = (double)(l+1)*(1.-eta2[l][1])*konst;
      double down = (double)l*(1.-eta2[l][0])*konst;
      SigmaAb[i-1] += down;
      SigmaAb[i] = up;
      //cout << l << " " << (up+down)/konst/(double)(2*l+1) << endl;
      if (parity == 1)
	{
	  SigmaAb_pos[i-1] += down;
	  SigmaAb_pos[i] = up;
	  SigmaAb_neg[i] = 0.;
	}
      else
	{
	  SigmaAb_neg[i-1] += down;
	  SigmaAb_neg[i] = up;
	  SigmaAb_pos[i] = 0.;
	}
      tot += up + down;
    }
  return tot;
  
}
//*************************************************************
  /**
   * return the total Elastic cross section, 
   * infinite for protons (don't try), finite for neutrons in mb
   * The function integrateWave() must be called first
   */

double scat::ElasticXsection()
{
  complex<double> one(1.,0.);
  //l==0 contribution
  double tot = pow(abs(one-eta[0][1]),2);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*pow(abs(one-eta[i][1]),2) + 
	(double)i*pow(abs(one-eta[i][0]),2);
    }
  return tot*konst;
  
} 


//*************************************************************
  /**
   * returns the total cross section, infinite for protons (don't try),
   * finite for neutrons in mb
   * The function integrateWave() must be called first
   */
double scat::TotXsection()
{
  //l==0 contribution
  double tot = 1+real(eta[0][1]);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*(1.-real(eta[i][1])) +
      (double)i*(1.-real(eta[i][0]));
    }
  return tot*konst*2.;
  
} 

//**************************************************************
  /**
   * returns the Rutherford differential cross section in mb
\param theta is the center-of-mass scattering angle in radians
  */
double scat::Rutherford(double theta)
{
  return 10.* pow(gamma,2)/4./Kwave2/pow(sin(theta/2.),4);
}
//*********************************************************************
  /**
   * After integration the wavefunction to rStop, this function returns 
   *Wavefunction magnitude and its derivative
   \param l0 is the orbital angular momentum of the nucleon
   \param j0 is the total angular momentum of the nucleon
   */
valarray<double> scat:: IntegrateBoundM( double j0, int l0)
{
 
  l = l0;
  j = j0;

  // initialize the coupled differential equation solver
  initMerson(.0001,.00001,.1,0);

  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  SpinOrbit->setAM(l,j);
  valarray<double> WaveFunct(2);

  // find initial and matching wave functions
  //initialise wave function
  // potential at start
  double Vstart = RealPotential(rStart);

  //derivative of potential at start
  double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;

  // initialize wavefunctions
  double fact = pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[0] = pow(rStart,l+1) - muhbar*(energyCM-Vstart)*fact; 
  // derivative of wavefunction
  fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
  WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                      - muhbar*(energyCM-Vstart)*fact;
  fact = pow(rStart,l+3)*muhbar/2./(double)(2*l+3);
  WaveFunct[1] += dVstart*fact;

  //integrate wavefunctions out to matching radius
  solveMerson(&WaveFunct,rStart,rStop);

  if (ok == 0) 
     {
       cout << "j= " << j << " l= " << l << " Ecm= " << energyCM << endl;
       valarray<double> finish(0);
       return finish;
     }

  return WaveFunct;
}
//*****************************************************************
  /**
   *calculates the Wavefunction and stores it in an array WaveArray
   */

void scat::GetWaveFunctionArrayM(double j0, int l0,
     double& derivative)
{

  j = j0;
  l = l0;
  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  SpinOrbit->setAM(l,j);

  WaveArray[0] = 0.;

 // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);
  valarray<double> WaveFunct(2);



  //initialise wave function
 // potential at start
  double Vstart = RealPotential(rStart);
     
  //derivative of potential at start
  double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;

  // initialize wavefunctions
  double fact = pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[0] = pow(rStart,l+1) - muhbar*(energyCM-Vstart)*fact; 
  // derivative of wavefunction
  fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
  WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                      - muhbar*(energyCM-Vstart)*fact;
  fact = pow(rStart,l+3)*muhbar/2./(double)(2*l+3);
  WaveFunct[1] += dVstart*fact;


  //integrate wavefunctions out to matching radius, but storing value
  // at intervales of deltaR

  double r1 = rStart;
  double r2 = r1+deltaR;
  WaveArray[0] = WaveFunct[0];


  for (int i=1;i<mWave;i++)
    {
     solveMerson(&WaveFunct,r1,r2);
     if (WaveFunct.size() == 0)
       {
	 cout << "problem in GetWaveFunctionArray " << endl;
	 abort();
       }
     WaveArray[i] = WaveFunct[0];
     WaveFunct = WaveFunct;



     r1 = r2;
     r2 += deltaR;
    }

  derivative = WaveFunct[1];
}
//*************************************************************
  /**
   * Normalizes the wavefunction
   */

void scat::normalizeWaveFunction(double E, double Efermi)
{

      double SumBar = 0.;
      double r = rStart;
      for(int kk=0;kk<nWave;kk++)
	{
	  WaveBar[kk] = WaveArray[kk];
          double delta = deltaR;
          if (kk == 0 || kk == nWave-1) delta /=2.;
	  SumBar += pow(WaveBar[kk],2)*delta;
          r += deltaR;
	}

}



//*********************************************************************
  /**
   * initialize integration of wavefunction
   */

void scat::initIntegration()
{
  rStart = .05;
  rStop = 12.;
  mWave = Nnumerov+1;
  nWave = mWave + 120;
  lMax = 60;
  llMax = lMax + 1;
// steps for wavefunction array in fm
  deltaR = (rStop-rStart)/(double)(mWave-1); 
  WaveArray = new double[nWave];
  WaveBar = new double[nWave];
  WaveMom = new double[100];

  Sigma = new double[llMax];
  eta2 = new double * [llMax];
  eta = new complex<double> * [llMax];
  phaseShift = new double * [llMax];
  for (int i=0;i<llMax;i++)
    {
      eta2[i] = new double [2];
      eta[i] = new complex<double> [2];
      phaseShift[i] = new double[2];
    }
 
  //set up array of Log derivatives.
  Nl = 7;
  LogDerMin = -123;
  LogDerMax = (int)EminRel;
  NLogDer = LogDerMax - LogDerMin + 1;
  
  LogDerArray = new double *[Nl+1];
  for (int i=0;i<=Nl;i++) LogDerArray[i] = new double [NLogDer];


  cout << " Making LogDer array" << endl;
  for (l = 0;l<=Nl;l++)
    {

     for (int i=LogDerMin;i<=LogDerMax;i++)
      {
        energyCM = (double) i;
        Kwave2 = kconstant*mu*energyCM;
        muhbar = kconstant*mu;
        Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1
        gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter

        LogDerArray[l][i-LogDerMin] = exteriorLogDer(l);
        
      }
    }
  cout << " Finished LogDer array" << endl;

  //make array for absorption cross section
  SigmaAb = new double[lMax+1];
  SigmaAb_pos = new double [lMax+1];
  SigmaAb_neg = new double [lMax+1];
}
//*************************************************************************
  /**
   *Destructor
   */
scat::~scat()
{

  if (flag == 0) return;
  delete [] WaveArray;
  delete [] WaveBar;
  delete [] WaveMom;
  delete [] Sigma;

  for (int i=0;i<llMax;i++)
    {
      delete [] eta2[i];
      delete [] eta[i];
      delete [] phaseShift[i];
    }
  delete [] eta2;
  delete [] eta;
  delete [] phaseShift;

  for (int i=0;i<=Nl;i++) delete [] LogDerArray[i];

  delete [] LogDerArray;

  delete [] SigmaAb;
  delete [] SigmaAb_pos;
  delete [] SigmaAb_neg;
}
//************************************************************************
  /**
   * calculates the log derivative of the exterior wavefunction
   \param l is the orbital angular momentum of the nucleon
  */
double scat::exteriorLogDer(int l)
{
  if (proton == 0)
    {
      sphericalB sph;
      if (energyCM < 0.) return sph.LogDer_k(l,rStop*Kwave)*Kwave;
      else if (energyCM >0.) return sph.LogDer_y(l,rStop*Kwave)*Kwave;
      else return -1.e32;
    }
  else
    {

      if (energyCM <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  double outwave,DoutwaveDr;
	  if (energyCM < 0.)
	    {
	      outwave = whittaker.getWaveFunction(gamma,l,rStop*Kwave);
	      DoutwaveDr = whittaker.derivative*Kwave;
	    }
	  else 
	    {
	      if (proton == 0) return 1.e32;
	      outwave = whittaker.ZeroEnergy( 2.*mu*Z*Zp/28.820*rStop,l);
	      DoutwaveDr = whittaker.derivative*2.*mu*Z*Zp/28.820;
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
//************************************************************************
  /**
   * adds the exterior part of the wavefunction to the array WaveArray
   * containing the interior part.
   \param l is the orbital angular momentum of the nucleon
  */
void scat::exteriorWaveFunct(int l)
{
  ANC=0.;
  if (proton == 0)
    {
      sphericalB sph;
      for (int i=0;i<=nWave-mWave;i++)
	{
	 double r =rStop + (double)i*deltaR;
	 double outwave;
         if (energyCM <= 0.) outwave = sph.k(l,r*Kwave);
	 else outwave = sph.y(l,r*Kwave);
	 if (i==0) ANC = WaveArray[mWave-1]/outwave;
	 else WaveArray[mWave+i-1] = outwave*ANC;
	}
    }
  else
    {
      if (energyCM <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*deltaR;
	      double outwave;
	      if (energyCM < 0.) outwave = 
              whittaker.getWaveFunction(gamma,l,r*Kwave);

	      else outwave = whittaker.ZeroEnergy( 2.*mu*Z*Zp/28.820*r,l);
	      if( i==0) ANC = WaveArray[mWave-1]/outwave;
	      else WaveArray[mWave+i-1] = outwave*ANC;
	    }
	}
      else  // from irregular Coulomb wave function for quasi-bound states
	{
	  coul Coulomb;
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*deltaR;
	      Coulomb.init(l,gamma,r*Kwave);
	      if( i==0) ANC = WaveArray[mWave-1]/Coulomb.G;
	      else WaveArray[mWave+i-1] = Coulomb.G*ANC;
	    }

	}
    }
}
//************************************************************
void scat::list()
{

  for (int i=-10;i<10;i++)
    {
      if (i==0) continue;
      energyCM = (double) i/4.;
      Kwave2 = kconstant*mu*energyCM;
      muhbar = kconstant*mu;
      Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1
      gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter

    }
}
//****************************************************************
double scat::LogDerDifference( double j, int l)
{
  //this would be fine if you didn't need the wavefunction
  //valarray<double> Waveff = IntegrateBound(j,l);


  //if wavefunction is needed must use this subroutine
  double Waveff[2];
  GetWaveFunctionArray(j, l,Waveff[1]);
  Waveff[0] = WaveArray[mWave-1];

  //log derivative of interior wave function
  //double LOgDerInterior = Waveff[1]/Waveff[0];

  //interpolate the exterior LogDerative
  double en = floor(energyCM);
  double delta = energyCM-en;
  int i = (int)en- LogDerMin;
  if (i < 0 || i > NLogDer) cout << " outside of LogDerArray, Ecm = " 
				<< energyCM << " i= " << i << endl;
  double LogDerExterior = LogDerArray[l][i];

  if (delta != 0.)
    {
      if (proton == 0 && fabs(energyCM)< 1.) 
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

  return Waveff[0] - Waveff[1]/LogDerExterior;
  
}
//**********************************************************************
  /**
   *returns the expectation value of the imaginary part of the mean field
   *used for spectral functions
   */
double scat::spectralWidth()
{
  double sum = 0.;
  double r = rStart;
      for (int kk=0;kk<nWave;kk++)
	{
          double delta = deltaR;
	  if (kk == 0 || kk == nWave-1) delta /= 2.;
	  double pot = Volume->ImaginaryPot(r) 
          + Surface->ImaginaryPot(r)
          + SpinOrbit->ImaginaryPotential(r);
	  sum += pot*pow(WaveBar[kk],2)*delta; 
          r += deltaR;
	}
  return sum;
}
//**************************************************************
  /**
   *returns the transmission coefficients
   *the function integrateWave() must be executed first
   \param l is the orbital angular momentum of the nucleon
   \param j is the total angular momentum of the nucleon
   */
double scat::TransCoef(int l, double j)
{
  if (l > lStop) return 0.;
  int jj = (int)j;
  if (jj == l)return (1.-eta2[l][1]);
  else return 1. - eta2[l][0];
}
//********************************************************************
//********************************************************************
  /**
   * calculates the integrated potentials and the RMS values of the 
   * potentials
   */

void scat::VIntegrals()
{
  double deltaR = .1;
  double sumReal = 0.;
  double sumImag = 0.;
  double sumSO = 0.;
  double sumRealR2 = 0.;
  double sumImagR2 = 0.;
  double sumSOR2 = 0.;
  SpinOrbit->setAM(0,0.5);
  for (int i = 0;i<200;i++)
    {
      double r = ((double)i+0.5)*deltaR;
      double potReal0 = HartreeFock->potential0(r);
      double potReal1 = HartreeFock->potential1(r);
      double potRealDis =
      + Volume->DispersiveCorrection(r)
	+ Surface->DispersiveCorrection(r);

      double potImag = Volume->ImaginaryPot(r) 
        + Surface->ImaginaryPot(r);


      double potSO = SpinOrbit->RealPotential(r); 


      if (beta_nl_R0 == 0)
      sumReal += potReal0*4.*pi*pow(r,2)*deltaR;
      else sumReal += potReal0*sqrt(pi)*beta_nl_R0/6.
            *(47./2.*pow(r,3)*exp(-pow(r/beta_nl_R0,2))
	    + 3.*r*pow(beta_nl_R0,2)*exp(-pow(r/beta_nl_R0,2)/2.)
            - 3.*beta_nl_R0*(-8.*pow(r,2)+pow(beta_nl_R0,2))
			      *sqrt(pi)*erf(r/2./beta_nl_R0))*deltaR;
      //see page 38 of my nonlocal notes

      if (beta_nl_R1 == 0)
      sumReal += potReal1*4.*pi*pow(r,2)*deltaR;
      else sumReal += potReal1*sqrt(pi)/beta_nl_R1
            *(47./2.*pow(r,3)*exp(-pow(r/beta_nl_R1,2))
	    + 3.*r*pow(beta_nl_R1,2)*exp(-pow(r/beta_nl_R1,2)/2.)
            - 3.*beta_nl_R1*(-8.*pow(r,2)+pow(beta_nl_R1,2))
			      *sqrt(pi)*erf(r/2./beta_nl_R1))*deltaR;
      //see page 38 of my nonlocal notes



      if (beta_nl_I == 0)
	{
         sumImag += potImag*4.*pi*pow(r,2)*deltaR;
	 sumReal += potRealDis*4.*pi*pow(r,2)*deltaR;
	}
      else
	{
	  double dv = sqrt(pi)/beta_nl_I/6.
             *(47./2.*pow(r,3)*exp(-pow(r/beta_nl_I,2))
	    + 3.*r*pow(beta_nl_I,2)*exp(-pow(r/beta_nl_I,2)/2.)
            - 3.*beta_nl_I*(-8.*pow(r,2)+pow(beta_nl_I,2))
			      *sqrt(pi)*erf(r/2./beta_nl_I))*deltaR;
         sumImag += potImag*dv;
	 sumReal += potRealDis*dv;
	}
      sumSO   += potSO*4.*pi*pow(r,2)*deltaR;

      sumRealR2 += potReal0*4.*pi*pow(r,4)*deltaR;
      sumImagR2 += potImag*4.*pi*pow(r,4)*deltaR;
      sumSOR2 += potSO*4.*pi*pow(r,4)*deltaR;

    }


  JReal = sumReal/A;
  JImag = sumImag/A;
  JSO = sumSO/pow(A,1./3.);

  RrmsReal = sqrt(sumRealR2/sumReal);
  RrmsImag = sqrt(sumImagR2/sumImag);  
  RrmsSO = sqrt(sumSOR2/sumSO);  
}
//******************************************************************
  /**
   * calculates the wavefunction at the specified momentum
   \param momentum is the momentum
   */
double scat::MomWave(double momentum)
{
  double r = rStart;
  double sum = 0.;
  for (int i = 0;i<nWave;i++)
    {
      double omega = momentum*r*.0050677;
      sum += WaveBar[i]*sin(omega)*deltaR;
      r += deltaR;
    }
  return 4./sqrt(pi)/momentum*sum;
}
//****************************************************************
  /**
   * calcuates the wavefunction in momentum space
   */

void scat::GetMomWave()
{
  double deltaMom = 3.;
  double sum = 0.;
  for (int i=0;i<100;i++)
    {
      double momentum = (double)i*deltaMom+.5;
      WaveMom[i] = MomWave(momentum);
      sum += WaveMom[i]*pow(momentum,2)*deltaMom;
    }
  sum *= 4.*pi;
  for (int i=0;i<100;i++)
    {
     WaveMom[i] /= sum;
     //cout << (double)i*deltaMom + .5 << " " << pow(WaveMom[i],2) << endl;
    }

}
//****************************************************
  /**
   *determines the number of nodes in the wavefunction
   */
int scat::Nodes()
{
  int nodes = 0;
  for (int i = 1;i<mWave-2;i++)
    {
 
      //cout << i << " " << WaveArray[i] << endl;
      if (WaveArray[i]*WaveArray[i-1] <= 0.) nodes++;

    }
  return nodes;
}
//******************************************************
//*****************************************************************
  /**
   *returns the absorption xsec associated with the isovector surface 
   *imaginary potential
   */
double scat::isovectorAbsorb()
{


  absorbAll = 0.;  //absorption xsec for all imaginary component
  absorbSurface = 0.;  //absorption xsec for surface imaginary component
  absorbVolume = 0.;
  absorbStandard = 0.;
  absorbSO = 0.; //spinOrbit

  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave; 


  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];

  valarray<double> WaveFunct(4);
  double Wave2Array[mWave];  //array containing the modulus squared of wavefunction

  l = 0;
  lStop = -1;
  for (;;) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;
          LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
	  SpinOrbit->setAM(l,j);

          // potential at start
          double Vstart = RealPotential(rStart);
          double Wstart = ImaginaryPotential(rStart);
          //derivative of potential at start
          double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;
          double dWstart = (ImaginaryPotential(rStart+.01) - Wstart)/0.01;

	  // initialize wavefunctions
          double fact = pow(rStart,l+3)/2./(double)(2*l+3);

          WaveFunct[0] = pow(rStart,l+1) 
                                - muhbar*(energyCM-Vstart)*fact; // real part
          WaveFunct[2] =  Wstart*muhbar*fact;              //imaginary part

          // derivative of wavefunction
          fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
          WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                           - muhbar*(energyCM-Vstart)*fact; // real
          WaveFunct[3] = muhbar*Wstart*fact;          // imaginary

          fact = muhbar*pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[1] += dVstart*fact;
          WaveFunct[3] += dWstart*fact;


          Wave2Array[0] = pow(WaveFunct[0],2)+pow(WaveFunct[2],2);          
          //integrate wavefunctions out to matching radius, but storing value
          // at intervales of deltaR
          double r1 = rStart;
          double r2 = r1+deltaR;
          for (int i=1;i<mWave;i++)
             {
              solveMerson(&WaveFunct,r1,r2);
              if (WaveFunct.size() == 0)
                {
	          cout << "problem in GetWaveFunctionArray " << endl;
	          abort();
                }
              Wave2Array[i] = pow(WaveFunct[0],2)+pow(WaveFunct[2],2);
              WaveFunct = WaveFunct;
              r1 = r2;
              r2 += deltaR;
             }

          if (ok == 0) 
	    {
              cout << "j= " << j << " l= " << l << " Ecm= " <<
		     energyCM << endl;
              return 0;
	    }

          //outWave gives derivates with respect to rho = Kwave*r
	  //but at this point WaveFunctOut have derivative with respect to r
          // so make them with respect to rho=kwave*r

	  WaveFunct[1] /= Kwave;
          WaveFunct[3] /= Kwave;


	  // match wave functions 
	  //real WaveFunct = AA*F + BB*G
	  double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
	  double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

	  // imaginary part => Wavefunct  = CC*F + DD*G
	  double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
	  double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];


	  double denominator = pow(AA+DD,2) + pow(CC-BB,2);
	  double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double etaImag = 2.*(AA*BB+CC*DD)/denominator;
          double eta2 = pow(etaReal,2) + pow(etaImag,2);
         

          //normalization factor
	  //complex<double> norm ((BB-CC)/2.,(AA+DD)/2.); 
          double norm = pow(BB-CC,2)+pow(AA+DD,2);

          double absorb_all = 0.;
          double absorb_S = 0.;
          double absorb_V = 0.;
          double absorb_SO = 0.;
	  double r = rStart;
          for (int i=0;i<mWave;i++)
	    {
              absorb_all += Wave2Array[i]/norm*deltaR*
                     ImaginaryPotential(r);
              absorb_S += Wave2Array[i]/norm*deltaR*
                      Surface->ImaginaryPot(r)*gammaRel;
              absorb_V += Wave2Array[i]/norm*deltaR*
		(Volume->ImaginaryPot(r))
                *gammaRel;
              if (r> 0)absorb_SO += Wave2Array[i]/norm*deltaR*
		SpinOrbit->ImaginaryPotential(r)*gammaRel;


	      /*
	      if (energyLab == 20. && l == 4 && plusMinus ==1)
		cout << r << " " << RealPotential(r) << endl;
	       
		cout << r << " " << Surface->ImaginaryPot(r) << 
		  " " << Wave2Array[i]/norm << endl;
		*/

	      r+= deltaR;


	    }

	  if (plusMinus == 1)
	    {
	      absorbAll += absorb_all*(double)(l+1);
	      absorbSurface += absorb_S*(double)(l+1);
	      absorbVolume += absorb_V*(double)(l+1);
              absorbSO += absorb_SO*(double)(l+1);
              absorbStandard += (1.-eta2)*konst*(double)(l+1);

              //if (energyLab == 20.) cout << 1 << " " << l << " "
              // << absorb_S*(double)(l+1) << endl; 

	    }
	  else 
	    {
	      absorbAll += absorb_all*(double)(l);
	      absorbSurface += absorb_S*(double)(l);
	      absorbVolume += absorb_V*(double)(l);
              absorbSO += absorb_SO*(double)(l);
              absorbStandard += (1.-eta2)*konst*(double)(l);

	     
              //if (energyLab == 20.) cout << 2 << " " << l << " "
              //    << absorb_S*(double)(l+1) << endl; 
	      

	    }

	}

      l++;
      if (l > lMax) break;
    }
  //the following are true classically
  //absorbIsovector *= -6.046*mu/pow(Kwave,3);
  //absorbIsoscaler *= -6.046*mu/pow(Kwave,3);
  //absorbVolume *= -6.046*mu/pow(Kwave,3);
  //absorbAll *= -6.046*mu/pow(Kwave,3);
  absorbSurface *= -4.*pi/Kwave/energyCM*10;
  absorbVolume *= -4.*pi/Kwave/energyCM*10;
  absorbAll *= -4.*pi/Kwave/energyCM*10;
  absorbSO *= -4.*pi/Kwave/energyCM*10;
 
   
  return absorbSurface;
}
  //*********************************************************************

  double scat::getSmatrixEikonal(double b)
    {
      LdotSigma = 0.;
      SpinOrbit->setAM(0,0.5);
      double const deltaZ = .1;
      double z = 0.;
      double maxR = 0.;
      double maxI = 0.;
      double sumR = 0.;
      double sumI = 0.;
      for (;;)
        {
          double r = sqrt(pow(b,2)+pow(z,2)); 
          double addR = HartreeFock->potential(r)
	  + Volume->DispersiveCorrection(r)
	    + Surface->DispersiveCorrection(r);

          maxR = max(maxR,abs(addR));
          if (z == 0.) addR /= 2;
          sumR += addR*deltaZ;

          double addI = ImaginaryPotential(r)/gammaRel;
          maxI = max(maxI,abs(addI));
          if (z == 0.) addI /= 2;
          sumI += addI*deltaZ;  

          if (abs(addR) < maxR/100. && abs(addI) < maxI/100.) break;
          if (z > 30.) break;
          z+= deltaZ;
        }

      double x = (pow(energyLab/m0+1,2)-1.);
      double beta = sqrt(x/(1.+x));
     


     double hvelocity = beta*197.326;
     double chiReal = -2.*sumR/hvelocity;
     double chiImag = -2.*sumI/hvelocity;


     //add the Coulomb term
     if (Zp == 1) chiReal += 2.*Z*e2/hvelocity*log(Kwave*b);


     //double phase = chiReal;
     double mag = exp(-chiImag);

     Sreal = mag*cos(chiReal);
     Simag = mag*sin(chiReal);


     return 1. - pow(mag,2);
    }

//*********************************************************************
//*****************************************************************
  /**
   * For positive energies, this Integrates the wave functions out from 
   * the origin to the matching radius, where the wavefunction is match to
   * Coulomb wavefunctions or spherical Bessel functions.
   * It determines the scattering phase shift. Uses the Numerov method for 
   * the integration
   */

int scat::integrateWaveN()
{


  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave;
  initNumerov(rStart,rStop,Nnumerov);

  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];


  l = 0;
  lStop = -1;
  for (;;) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;
          LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
	  SpinOrbit->setAM(l,j);



          nonLocalIntegration(complex<double>(pow(rStart,l+1),0.),
		       complex<double>((double)(l+1)*pow(rStart,l),0.));

	  

	  complex<double> y_end = y[Nnumerov];
          complex<double> dydr_end = (y[Nnumerov+1] - y[Nnumerov-1])/2./dx; 


	  //rescale in case the numbers are too big, 
	  //this happend for large l with large beta sometimes.
	  //rescaling stops overflow latter in the code
          dydr_end /= y_end.real();
          y_end /= y_end.real();




          //outWave gives derivates with respect to rho = Kwave*r
	  //but at this point WaveFunctOut have derivative with respect to r
          // so make them with respect to rho=kwave*r
	  dydr_end /= Kwave;




	  // match wave functions 
	  //real WaveFunct = AA*F + BB*G
	  double  BB = outWave.dF[l]*y_end.real() 
	    - outWave.F[l]*dydr_end.real();
	  double AA = -outWave.dG[l]*y_end.real() 
	    + outWave.G[l]*dydr_end.real();

	  // imaginary part => Wavefunct  = CC*F + DD*G
	  double DD = outWave.dF[l]*y_end.imag() 
	    - outWave.F[l]*dydr_end.imag();
	  double CC = -outWave.dG[l]*y_end.imag() 
	    + outWave.G[l]*dydr_end.imag();


	  double denominator = pow(AA+DD,2) + pow(CC-BB,2);
	  double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double etaImag = 2.*(AA*BB+CC*DD)/denominator;

          int up = (plusMinus+1)/2;  //=0 form down and =1 for up
	  phaseShift[l][up] = atan2(etaImag,etaReal);
          if (phaseShift[l][up] < 0.) phaseShift[l][up] += 2.*pi;
          phaseShift[l][up] /=2.;
         
	  eta2[l][up] = pow(etaReal,2) + pow(etaImag,2);
          eta[l][up] = complex<double>(etaReal,etaImag);





	}
      //check to see if we have included enough l-waves
      if (1.-eta2[l][1] < 0.001 && 1.-eta2[l][0] < 0.001 )
	{
          lStop = l;
	  break;
	}
      if (isnan(eta2[l][1]) || isnan(eta2[l][0]))
	 {
	   lStop = l-1;
           //cout << 1. - eta2[lStop][1] << " " << 1. - eta2[lStop][0] << endl;
	   break;
	 } 
      l++;
      if (l > lMax) break;
    }

      if (lStop == -1) cout << "increase lMax" << endl;
      return 1;
}
//***************************************************************
int scat::integrateWave()
{
  if (Imethod == 0) return integrateWaveM();
  else return integrateWaveN();
}
//****************************************************************
double scat::nonlocalFact(double r1, double r2, double beta)
{

  if ( l == 0)
    return (exp(-pow((r1-r2)/beta,2)) - exp(-pow((r1+r2)/beta,2)))/sqrt(pi)/beta;
  double x = 2.*r1*r2/pow(beta,2);
  if ( l == 1) 
    return (exp(-pow((r1-r2)/beta,2))*(1.-1./x) + 
          exp(-pow((r1+r2)/beta,2))*(1.+1./x))/sqrt(pi)/beta;
  if ( l == 2)
    return (exp(-pow((r1-r2)/beta,2))*(1.-3./x + 3./pow(x,2)) + 
	    exp(-pow((r1+r2)/beta,2))*(-1.-3./x-3./pow(x,2)))/sqrt(pi)/beta;   
  if ( l == 3)
    return (exp(-pow((r1-r2)/beta,2))*(1.-6./x + 15./pow(x,2)-15./pow(x,3)) + 
	    exp(-pow((r1+r2)/beta,2))*(1.+6./x + 15./pow(x,2)+15./pow(x,3)))
            /sqrt(pi)/beta;   

  if ( l == 4)
    return (exp(-pow((r1-r2)/beta,2))*(1.-10./x + 45./pow(x,2)-105./pow(x,3)
                                       +105./pow(x,4)) + 
	    exp(-pow((r1+r2)/beta,2))*(-1.-10./x-45./pow(x,2)-105./pow(x,3)-
				       105./pow(x,4)))
            /sqrt(pi)/beta;   

  if ( l == 5)
    return (exp(-pow((r1-r2)/beta,2))*(1.-15./x + 105./pow(x,2)-420./pow(x,3)
                                       +945./pow(x,4)-945./pow(x,5)) + 
	    exp(-pow((r1+r2)/beta,2))*(1.+15./x+105./pow(x,2)+420./pow(x,3)+
				       945./pow(x,4)+945./pow(x,5)))
            /sqrt(pi)/beta;   


  if ( l == 6)
    return (exp(-pow((r1-r2)/beta,2))*(1.-21./x + 210./pow(x,2)-1260./pow(x,3)
                            +4725./pow(x,4)-10395./pow(x,5)+10395./pow(x,6)) + 
	    exp(-pow((r1+r2)/beta,2))*(-1.-21./x-210./pow(x,2)-1260./pow(x,3)-
			     4725./pow(x,4)-10395./pow(x,5)-10395./pow(x,6)))
            /sqrt(pi)/beta;
  if ( l == 7)
    return (exp(-pow((r1-r2)/beta,2))*(1.-28./x + 378./pow(x,2)-3150./pow(x,3)
                            +17325./pow(x,4)-62370./pow(x,5)+135135./pow(x,6)-
				       135135./pow(x,7)) + 
	    exp(-pow((r1+r2)/beta,2))*(1.+28./x+378./pow(x,2)+3150./pow(x,3)+
			     17325./pow(x,4)+62370./pow(x,5)+135135./pow(x,6)
				       +135135./pow(x,7)))
            /sqrt(pi)/beta;

  if ( l == 8)
    return (exp(-pow((r1-r2)/beta,2))*(1.-36./x + 630./pow(x,2)-6930./pow(x,3)
                            +51975./pow(x,4)-270270./pow(x,5)+945945./pow(x,6)-
				       2027025./pow(x,7)+ 2027025./pow(x,8)) + 
	    exp(-pow((r1+r2)/beta,2))*(-1.-36./x-630./pow(x,2)-6930./pow(x,3)-
			     51975./pow(x,4)-270270./pow(x,5)-945945./pow(x,6)
				       -2027025./pow(x,7)-2027025./pow(x,8)))
            /sqrt(pi)/beta;


  if ( l == 9)
    return (exp(-pow((r1-r2)/beta,2))*(1.-45./x + 990./pow(x,2)-13860./pow(x,3)
                  +135135./pow(x,4)-945945./pow(x,5)+4729725./pow(x,6)
		  -16216200./pow(x,7)+ 34459425./pow(x,8)-34459425/pow(x,9)) + 
	    exp(-pow((r1+r2)/beta,2))*(1.+45./x+990./pow(x,2)+13860./pow(x,3)
		  + 135135./pow(x,4)+945945./pow(x,5)+4729725./pow(x,6)
		  + 16216200./pow(x,7)+34459425./pow(x,8)+34459425./pow(x,9)))
            /sqrt(pi)/beta;


  sphericalB sph;

  if (x < 400.)
    {

    sph.I(l,x);  //modified spherical bessels function of the first kind
    double out = 4.*r1*r2/sqrt(pi)/pow(beta,3)*
    exp(-(pow(r1,2)+pow(r2,2))/pow(beta,2));
    return out*sph.II[l];
    }
  else 
    {
    sph.asymptoticI(l,x);  //calculates the modified 
                           //spherical bessels function of the first kind
                           //times the factor exp(-x)/2./x 
                           //in the asymptotic limit 
    double out = sph.II[l]/sqrt(pi)/beta*
                 exp(-pow((r1-r2)/beta,2));

    if (isnan(out)|| isinf(out))
      {
        cout << "nonLocalFactor problem, r1= "<< r1 << " r2= " << r2 << 
	  "beta = " << beta << endl;
      }
   return out;
    }

}
//*****************************************************************
void scat::nonLocalArray1()
{

  //source term and local potentials

  for (int i=0;i<=Nnumerov+1;i++)
    {

      //get F factor
      double F = 0;
      if (i > 0)
	{
         complex<double> slope;

         //estimate slope of wavefunction   
         if (i < Nnumerov+1) slope = (y[i+1]- y[i-1])/2./dx;
         else slope = (3.*y[i] + y[i-2] - 4.*y[i-1])/2./dx;


         complex<double> stuff= y[i]/slope;
      
         F = exp(-100.*pow(abs(stuff),2));
         stuff = pow(x[i],l+1)/y[i] -1.;

         F *= 1. - exp(-100.*pow(abs(stuff),2));
	}

      // now do the convolution
      complex<double> pot(0.,0.);
      double r1 = x[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x[j];
          if (r1 - r2 > 3.*beta_max) continue;
          if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

	  if (beta_nl_R0 > 0. && ddr < beta_nl_R0)
	    pot += dx*nonlocalFact(r1,r2,beta_nl_R0)*y[j]*
	      complex<double>(HartreeFock->U0(r1,r2),0.);

	  if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0)
	    pot += dx*nonlocalFact(r1,r2,beta_nl_R1)*y[j]*
	      complex<double>(HartreeFock->U1(r1,r2),0.);



	  if (beta_nl_I > 0. && ddr < beta_nl_I)
	    {
	      complex<double> fs = Surface->U(r1,r2);
	      complex<double> fv = Volume->U(r2,r2);

	    pot += dx*nonlocalFact(r1,r2,beta_nl_I)*y[j]*(fv+fs);
	    }

	}


      source[i] = muhbar*F*pot*gammaRel;
      complex <double> vloc = muhbar*(1.-F)/y[i]*pot*gammaRel;



      // now add in the remainding local parts of the potential
      double fr = 0;
      if (proton == 1) // finally the Coulomb potential
        {
          if (x[i] > Rc) fr = e2*Z*Zp/x[i];
          else fr = e2*Z*Zp/2./Rc*(3.-pow(x[i]/Rc,2));
        }

      if (x[i] > 0.) 
        {
          // next the spin-orbit potential
         fr += SpinOrbit->RealPotential(x[i]);
        }

      fr += HartreeFock->potentialSurface(x[i]);

      if (beta_nl_R0 == 0.)
      fr += HartreeFock->potential0(x[i]);


      if (beta_nl_R1 == 0. && HartreeFock->A != 0.)
      fr += HartreeFock->potential1(x[i]);

      
      if (beta_nl_I == 0.)
	{
         double two = Volume->DispersiveCorrection(x[i]);
         double three = Surface->DispersiveCorrection(x[i]);
	 fr += two + three;
	}
      
       
      fr *= gammaRel*muhbar;
      if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);
      fr = -(Kwave2 - fr); 

      double fi = 0.;
       if (x[i] > 0.) 
         {
          // next the spin-orbit potential
           fi= SpinOrbit->ImaginaryPotential(x[i]);
          }

       if (beta_nl_I ==  0.) fi +=  Volume->ImaginaryPot(x[i])
	 + Surface->ImaginaryPot(x[i]);

       fi *= gammaRel*muhbar;


      f[i] = vloc +complex<double>(fr,fi);
 
    }


}
//*************************************************************
void scat::nonLocalIntegration(complex<double> y0, complex<double>dydr0)
  {
    nonLocalArray0();
    solveNumerov(y0,dydr0);
    if (local) return;
    if (beta_nl_R0 == 0. && beta_nl_R1 == 0. && beta_nl_I == 0) return;

    for (int itry = 1;itry<4;itry++)
      {
	nonLocalArray1();
        solveNumerov(y0,dydr0);

      }

    //if (l == 10)printWave();


  }
//*********************************************************
  /**
   * sets the source and F arrays for the the first
   * iteration of the wavefunction for the nonlocal calculation
   */
void scat::nonLocalArray0()
{

  double V = - HartreeFock->Vvol;
  for (int i=0;i<15;i++)
    {
     double  k2 = (energyCM - V - CoulDis)*kconstant*mu;
     V = -HartreeFock->Vvol/(1.+HartreeFock->A)*exp(-k2*pow(beta_nl_R0,2)/4.)
        - HartreeFock->A*HartreeFock->Vvol/(1.+HartreeFock->A)*
           exp(-k2*pow(beta_nl_R1,2)/4.) - CoulDis;

    }
  HartreeFock->setEquivLocalDepth(-V);

  double factorR0 = exp(-pow(beta_nl_R0*Kwave/2.,2));
  double factorR1 = exp(-pow(beta_nl_R1*Kwave/2.,2));
  double VV = (HartreeFock->Vvol*factorR0 
	       + HartreeFock->A*HartreeFock->Vvol*factorR1)/(1.+HartreeFock->A);
  HartreeFock->setEquivLocalDepth(VV);
  double factorI = exp(-pow(beta_nl_I*Kwave/2.,2));


    for (int i=0;i<N+2;i++)
      {
       double fr = HartreeFock->potentialEquivLocal(x[i]);

       double two = Volume->DispersiveCorrection(x[i])*factorI;
       double three = Surface->DispersiveCorrection(x[i])*factorI;
       fr += two + three;



       if (proton == 1) //the Coulomb potential
         {
            if (x[i] > Rc) fr += e2*Z*Zp/x[i];
            else fr += e2*Z*Zp/2./Rc*(3.-pow(x[i]/Rc,2));
         }

       if (x[i] > 0.) 
         {
          // next the spin-orbit potential
          fr += SpinOrbit->RealPotential(x[i]);
          }

       fr *= gammaRel*muhbar;

       if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);
       fr = -(Kwave2 - fr);

       double fi = (Volume->ImaginaryPot(x[i])+Surface->ImaginaryPot(x[i]))
                   *factorI;
       if (x[i] > 0.) 
         {
           // next the spin-orbit potential
           fi += SpinOrbit->ImaginaryPotential(x[i]);
         }
       fi *= gammaRel*muhbar;

       f[i] = complex<double>(fr,fi);

       source[i] = complex<double>(0.,0.);


      }

}

//********************************************************
void scat::printVloc()
{
  for (int i=0;i<=Nnumerov+1;i++)
    {
      // now do the convolution
      complex<double> pot(0.,0.);
      double r1 = x[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x[j];
          if (r1 - r2 > 3.*beta_max) continue;
          if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

	  if (beta_nl_R0 > 0.  && ddr < beta_nl_R0)
	    {
	      double fr = HartreeFock->U0(r1,r2);
	      pot += dx*nonlocalFact(r1,r2,beta_nl_R0)*y[j]*
	      complex<double>(fr,0.);
	    }

	  if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0.)
	    {
	      double fr = HartreeFock->U1(r1,r2);
	      pot += dx*nonlocalFact(r1,r2,beta_nl_R1)*y[j]*
	      complex<double>(fr,0.);
	    }


	  if (beta_nl_I > 0. && ddr < beta_nl_I)
	    {
	      complex<double> fs = Surface->U(r1,r2);
	      complex<double> fv = Volume->U(r1,r2);

	      pot += dx*nonlocalFact(r1,r2,beta_nl_I)*y[j]*(fs+fv);


	    }

	}
      complex<double> out = pot/y[i];
      cout << r1 << " " << out.real() << " " << out.imag()  << endl;
    }
  abort();
}
//***************************************************
void scat::printWave()
{
  for (int i=0;i<=Nnumerov+1;i++)
    {

      // now do the convolution
      complex<double> pot(0.,0.);
      double r1 = x[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x[j];
          if (r1 - r2 > 3.*beta_max) continue;
          if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

	  if (beta_nl_R0 > 0. && ddr < beta_nl_R0)
	    {
	      double fr = HartreeFock->U0(r1,r2);
	      pot += dx*nonlocalFact(r1,r2,beta_nl_R0)*y[j]*
	      complex<double>(fr,0.);
	    }

	  if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0.)
	    {
	      double fr = HartreeFock->U1(r1,r2);
	    pot += dx*nonlocalFact(r1,r2,beta_nl_R1)*y[j]*
	      complex<double>(fr,0.);
	    }
	  if (beta_nl_I > 0. && ddr < beta_nl_I)
	    {
	      complex<double> fs = Surface->U(r1,r2);
	      complex<double> fv = Volume->U(r1,r2);

	      pot += dx*nonlocalFact(r1,r2,beta_nl_I)*y[j]*(fs+fv);
	    }

	}


      // now add in the remainding local parts of the potential
      double fr = 0;
      if (proton == 1) // finally the Coulomb potential
        {
          if (x[i] > Rc) fr = e2*Z*Zp/x[i];
          else fr = e2*Z*Zp/2./Rc*(3.-pow(x[i]/Rc,2));
        }

      if (x[i] > 0.) 
        {
          // next the spin-orbit potential
         fr += SpinOrbit->RealPotential(x[i]);
        }

      fr += HartreeFock->potentialSurface(x[i]);

      if (beta_nl_R0 == 0.)
      fr += HartreeFock->potential0(x[i]);

      if (beta_nl_R1 == 0. && HartreeFock->A != 0.)
      fr += HartreeFock->potential1(x[i]);
  
      if (beta_nl_I == 0.)
      fr += Volume->DispersiveCorrection(x[i])
	+Surface->DispersiveCorrection(x[i]);


      fr *= gammaRel*muhbar;
      if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);
      fr = -(Kwave2 - fr); 

      double fi = 0.;
       if (x[i] > 0.) 
         {
          // next the spin-orbit potential
           fi= SpinOrbit->ImaginaryPotential(x[i]);
          }

       if (beta_nl_I ==  0.) fi +=  Volume->ImaginaryPot(x[i])
	 + Surface->ImaginaryPot(x[i]);

       fi *= gammaRel*muhbar;


      f[i] = pot*gammaRel*muhbar +complex<double>(fr,fi)*y[i];
      if (i == 0) continue;
      complex<double> d2 = (y[i+1]+y[i-1]-2.*y[i])/dx/dx;

      //cout << x[i] << " " << y[i].imag() << " " << f[i].imag() << " " << d2.imag() << endl;
      cout << x[i] << " " << y[i].real() << " " << y[i].imag() << endl;
    }
  abort();
   
}


//*********************************************************************
  /**
   * After integration the wavefunction to rStop, this function returns 
   *Wavefunction magnitude and its derivative
   \param l0 is the orbital angular momentum of the nucleon
   \param j0 is the total angular momentum of the nucleon
   */
valarray<double> scat:: IntegrateBoundN( double j0, int l0)
{
 
  l = l0;
  j = j0;

 

 // initialize the coupled differential equation solver
  initNumerovR(rStart,rStop,Nnumerov);

  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  SpinOrbit->setAM(l,j);

  double y0 = pow(rStart,l+1);
  double dydr0 = (double)(l+1)*pow(rStart,l);
  
  if (phase == 0)
    {
      nonLocalArray0R(); //use perey and Buck equiv local potential(1st iteration)
     solveNumerovR(y0,dydr0);
    }
  else if (phase == 1)
    {
      nonLocalArray1R(false); //use TELP equiv local potential(subsequent inerations)
     solveNumerovR(y0,dydr0);
    }
  else 
    nonLocalIntegration(y0,dydr0);


  valarray<double> out(2);
  out[0] =  y_R[Nnumerov];
  out[1] = (y_R[Nnumerov+1]-y_R[Nnumerov-1])/2./dx_R;
  return out;
}
//*************************************************************
void scat::nonLocalIntegration(double y0, double dydr0)
  {
    nonLocalArray0R();
    solveNumerovR(y0,dydr0);

    if (local) return;
    if (beta_nl_R0 == 0. && beta_nl_R1 == 0 && beta_nl_I == 0) return;


    for (int itry = 1;itry<8;itry++)
      {
	nonLocalArray1R(true);
        solveNumerovR(y0,dydr0);

      }
  }
//*********************************************************
  /**
   * sets the source and F arrays for the the first
   * iteration of the wavefunction for the nonlocal calculation
   */
void scat::nonLocalArray0R()
{




  double V = - HartreeFock->Vvol;
  for (int i=0;i<15;i++)
    {
     double  k2 = (energyCM - V - CoulDis)*kconstant*mu;
     V = -HartreeFock->Vvol/(1.+HartreeFock->A)*exp(-k2*pow(beta_nl_R0,2)/4.)
        - HartreeFock->A*HartreeFock->Vvol/(1.+HartreeFock->A)*
           exp(-k2*pow(beta_nl_R1,2)/4.) - CoulDis;

    }
  HartreeFock->setEquivLocalDepth(-V);

  double factorI = 1.;

    for (int i=0;i<N+2;i++)
      {
	double fr = HartreeFock->potentialEquivLocal(x_R[i]);

       double two = Volume->DispersiveCorrection(x_R[i])*factorI;
       double three = Surface->DispersiveCorrection(x_R[i])*factorI;
       fr += two + three;
       if (proton == 1) //the Coulomb potential
         {
            if (x_R[i] > Rc) fr += e2*Z*Zp/x_R[i];
            else fr += e2*Z*Zp/2./Rc*(3.-pow(x_R[i]/Rc,2));
         }

       if (x_R[i] > 0.) 
         {
          // next the spin-orbit potential
          fr += SpinOrbit->RealPotential(x_R[i]);
          }

       fr *= muhbar;

       if (l > 0) fr += (double)(l*(l+1))/pow(x_R[i],2);
       fr = -(Kwave2 - fr);


       f_R[i] = fr;

       source_R[i] = 0.;
      }
}
//*****************************************************************
void scat::nonLocalArray1R(bool generateVloc)
{

  //source term and local potentials
  if (generateVloc) generateLocEquiv();

  for (int i=0;i<=Nnumerov+1;i++)
    {

      // now add in the remainding local parts of the potential
      double fr = 0;
      if (proton == 1) // finally the Coulomb potential
        {
          if (x_R[i] > Rc) fr = e2*Z*Zp/x_R[i];
          else fr = e2*Z*Zp/2./Rc*(3.-pow(x_R[i]/Rc,2));
        }

      if (x_R[i] > 0.) 
        {
          // next the spin-orbit potential
         fr += SpinOrbit->RealPotential(x_R[i]);
        }

      fr += HartreeFock->potentialSurface(x_R[i]);

      if (beta_nl_R0 == 0.)
      fr += HartreeFock->potential0(x_R[i]);

      if (beta_nl_R1 == 0. && HartreeFock->A != 0.)
      fr += HartreeFock->potential1(x_R[i]);

      
      if (beta_nl_I == 0.)
	{
         double two = Volume->DispersiveCorrection(x_R[i]);
         double three = Surface->DispersiveCorrection(x_R[i]);
	 fr += two + three;
	}
      
       
      fr *= muhbar;
      if (l > 0) fr += (double)(l*(l+1))/pow(x_R[i],2);
      fr = -(Kwave2 - fr); 


      f_R[i] = vlocR[i] + fr;
 
    }


}
//*************************************************************

valarray<double> scat::IntegrateBound(double j, int l)
{
  if (Imethod == 0) return IntegrateBoundM( j, l);
  else return IntegrateBoundN(j, l);
}

//*****************************************************************
  /**
   *calculates the Wavefunction and stores it in an array WaveArray
   */

void scat::GetWaveFunctionArrayN(double j0, int l0,
     double& derivative)
{

  j = j0;
  l = l0;
  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  SpinOrbit->setAM(l,j);
  valarray<double> WaveFunct = IntegrateBoundN(j0,l0);


  for (int i=0;i<Nnumerov+1;i++)
    {
      WaveArray[i] = y_R[i];
    }

  derivative = WaveFunct[1];
}
//**************************************************
void scat::GetWaveFunctionArray(double j, int l, double & derivative)
{
  if (Imethod == 0) return GetWaveFunctionArrayM( j, l,derivative);
  else return GetWaveFunctionArrayN(j, l,derivative);
}
//********************************************************
void scat::printVlocR()
{
  cout << "Kwave2 = " << Kwave2 << endl;
  cout << "Kwave = " << Kwave << endl;
  cout << "Ecm = " << energyCM << endl;





  for (int i=0;i<=Nnumerov+1;i++)
    {
      // now do the convolution
      double pot = 0.;
      double r1 = x_R[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x_R[j];
          if (r1 - r2 > 3.*beta_max) continue;
          if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

	  if (beta_nl_R0 > 0. && ddr < beta_nl_R0)
	    {
	      double fr = HartreeFock->U0(r1,r2);
	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_R0)*y_R[j]*fr;

	    }

	  if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0.)
	    {
	      double fr = HartreeFock->U1(r1,r2);
	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_R1)*y_R[j]*fr;

	    }


	  if (beta_nl_I > 0. && ddr < beta_nl_I)
	    {
	      double two = (Volume->U(r1,r2)).real();
	      double three = (Surface->U(r1,r2)).real();
	      double fr = two + three;
	      pot += dx_R*nonlocalFact(r1,r2,beta_nl_I)*y_R[j]*fr;
	    }

	}
      double out = pot/y_R[i];
      cout << r1 << " " << y_R[i] << " " << out  << endl;
    }
  abort();
}
//***************************************************
void scat::printWaveR()
{
  cout << "Kwave2 = " << Kwave2 << endl;
  cout << "Kwave = " << Kwave << endl;
  cout << "Ecm = " << energyCM << endl;


  

  for (int i=0;i<=Nnumerov+1;i++)
    {

      // now do the convolution
      double pot=0.;
      double r1 = x_R[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x_R[j];
          if (r1 - r2 > 3.*beta_max) continue;
          if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

	  if (beta_nl_R0 > 0. && ddr < beta_nl_R0)
	    {
	      double fr = HartreeFock->U0(r1,r2);
	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_R0)*y_R[j]*fr;
	    }

	  if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0.)
	    {
	      double fr = HartreeFock->U1(r1,r2);
	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_R1)*y_R[j]*fr;
	    }
	  if (beta_nl_I > 0. && ddr < beta_nl_I)
	    {
	      double two = (Volume->U(r1,r2)).real();
	    double three = (Surface->U(r1,r2)).real();
            double fr = two + three;

	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_I)*y_R[j]*fr;
	    }

	}


      // now add in the remainding local parts of the potential
      double fr = 0;
      if (proton == 1) // finally the Coulomb potential
        {
          if (x_R[i] > Rc) fr = e2*Z*Zp/x_R[i];
          else fr = e2*Z*Zp/2./Rc*(3.-pow(x_R[i]/Rc,2));
        }

      if (x_R[i] > 0.) 
        {
          // next the spin-orbit potential
         fr += SpinOrbit->RealPotential(x_R[i]);
        }

      fr += HartreeFock->potentialSurface(x_R[i]);

      if (beta_nl_R0 == 0.)
      fr += HartreeFock->potential0(x_R[i]);

      if (beta_nl_R1 == 0. && HartreeFock->A != 0.)
      fr += HartreeFock->potential1(x_R[i]);
  
      if (beta_nl_I == 0.)
      fr += Volume->DispersiveCorrection(x_R[i])
	+Surface->DispersiveCorrection(x_R[i]);


      fr *= muhbar;
      if (l > 0) fr += (double)(l*(l+1))/pow(x_R[i],2);
      fr = -(Kwave2 - fr); 


      f_R[i] = pot*muhbar + fr*y_R[i];
      if (i == 0) continue;
      double d2 = (y_R[i+1]+y_R[i-1]-2.*y_R[i])/dx_R/dx_R;

      cout << x_R[i] << " " << f_R[i]  << " " << d2 << endl;

    }
  abort();
   
}
//*****************************************************************

void scat::generateLocEquiv()
{
  for (int i=0;i<=Nnumerov+1;i++)
    {
	
      //get F factor
      double F = 0;
      if (i > 0)
        {
          double slope;

          //estimate slope of wavefunction   
          if (i < Nnumerov+1) slope = (y_R[i+1]- y_R[i-1])/2./dx_R;
	  else slope = (3.*y_R[i] + y_R[i-2] - 4.*y_R[i-1])/2./dx_R;


          double stuff= y_R[i]/slope;
      
          F = exp(-100.*pow(stuff,2));
          stuff = pow(x_R[i],l+1)/y_R[i] -1.;

          F *= 1. - exp(-100.*pow(stuff,2));
        }

       // now do the convolution
      double pot = 0.;
      double r1 = x_R[i];
      for (int j=0;j<=Nnumerov+1;j++)
        {
	  double r2 = x_R[j];
	  if (r1 - r2 > 3.*beta_max) continue;
	  if (r2 - r1 > 3.*beta_max) break;
	  double ddr = abs(r1-r2)/3.;

          if (beta_nl_R0 > 0. && ddr < beta_nl_R0)
	    pot += dx_R*nonlocalFact(r1,r2,beta_nl_R0)*y_R[j]*
	      HartreeFock->U0(r1,r2);

	   if (beta_nl_R1 > 0. && ddr < beta_nl_R1 && HartreeFock->A != 0.)
	     pot += dx_R*nonlocalFact(r1,r2,beta_nl_R1)*y_R[j]*
	       HartreeFock->U1(r1,r2);

	   if (beta_nl_I > 0. && ddr < beta_nl_I)
	     {
	       double two =  (Volume->U(r1,r2)).real();
	       double three = (Surface->U(r1,r2)).real();
	       double fr = two + three;
	       pot += dx_R*nonlocalFact(r1,r2,beta_nl_I)*y_R[j]*fr;
	     }

	 }


       source[i] = muhbar*F*pot;
       vlocR[i] = muhbar*(1.-F)/y_R[i]*pot;
    }
}
