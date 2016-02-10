#include "scatterRspace.h"

double const scatterRspace::pi=acos(-1.);
double const scatterRspace::kconstant = .048192;
double const scatterRspace::m0 = 931.5; // nucleon mass
double const scatterRspace::e2 = 1.44; // Coulomb constant


/**
 * initialization for a new target and projectile
\param Z0 - target atomic number
\param Zp0 - projectile atomic number (0= neytron, 1=proton)
\param A0 - target atomic number
*/
void scatterRspace::init(double Z0, double Zp0, double A0)
{
  Z = Z0;
  A = A0;
  Zp = Zp0;
  mu = A/(1.+A);
}
//***************************************************************************
  /**
   * Constructor
   \param Z0 - target atomic number
   \param Zp0 - projectile atomic number (0= neytron, 1=proton)
   \param A0 - target atomic number
   \param Pot0 - pointer to object which calculates potential energy
   \param title0 - title for use in compound class (compound elastic contribution)
  */
scatterRspace::scatterRspace(double Z0, double Zp0, double A0, pot*Pot0, 
string *title0) :compound( title0)
{
  Pot = Pot0;
  init(Z0,Zp0,A0);
 // rStart = .05;
  rStart = .05;

  rStop = 10.;
  Nnumerov =400;
 // Nnumerov = 200;
  initNumerov(rStart,rStop,Nnumerov);

  lMax = 60;
 // lMax =200;

  llMax = lMax + 1;

  Sigma = new double[llMax];
  S2 = new double * [llMax];
  S = new complex<double> * [llMax];
  phaseShift = new double * [llMax];
  for (int i=0;i<llMax;i++)
    {
      S2[i] = new double [2];
      S[i] = new complex<double> [2];
      phaseShift[i] = new double[2];
    }

  //make array for absorption cross section
  SigmaAb = new double[lMax+1];
  SigmaAb_pos = new double [lMax+1];
  SigmaAb_neg = new double [lMax+1];

}

//************************************************************
  /**
   * destructor
   */
scatterRspace::~scatterRspace()
{
  delete [] Sigma;

  for (int i=0;i<llMax;i++)
    {
      delete [] S2[i];
      delete [] S[i];
      delete [] phaseShift[i];
    }
  delete [] S2;
  delete [] S;

  delete [] SigmaAb;
  delete [] SigmaAb_pos;
  delete [] SigmaAb_neg;
}
//*******************************************************************
  /**
   * Set the pointer to a new pot class (calculates the potential energy)
   \param Pot0 - calls that calculates the potential energy
  */
void scatterRspace::setPot(pot*Pot0)
{
  Pot = Pot0;
}
//****************************************************************
  /**
   * initialized the calculation for a new projectle energy
   \param Ecm0 - center of mass energy in MeV
   \param Elab0 - lab frame nergy in MeV
  */
void scatterRspace::newEnergy(double Ecm0, double Elab0)
{
  Ecm = Ecm0;
  Elab = Elab0;

  //relativistic verion
//  Kwave2 = kconstant*pow(A/(1.+ A + Ecm/m0),2)*Elab*(Elab/2./m0 + 1.);  
  // classical version
 // mu  = mu/ 1.0066;
 Kwave2 = kconstant*mu*Ecm;
  

 Kwave = sqrt(Kwave2);
  
// cout<<"Ecm =" << Ecm<<"\t Kwave in scatterRspace is" <<Kwave<<endl;
// cout<<"kconstant =" << kconstant << "mu is " << mu<<endl;
 //cout<< "Mass"<<kconstant*mu<<endl;
 // gammaRel = 2.*(Ecm/m0 +1.)/(Ecm/m0 + 2.);
  gammaRel = 1.;//;2.;//*(Ecm/m0 +1.)/(Ecm/m0 + 2.);

  muhbar = Kwave2/Ecm;
  //cout<<"muhbar "<<muhbar<<endl;

  gamma = fabs(gammaRel*Z*Zp*e2*Kwave/Ecm/2.); //Sommerfeld parameter

  konst = pi/pow(Kwave,2); //const for cross sections in units of Fermi squared
  konst *= 10.; //units of mb
}
//********************************************************************
  /**
   * Solves the Schrodinger Equation for each l and J value and 
   * calculates the S matrix
   \param Ecm0 - center of mass energy in MeV
   \param Elab0 - lab frame nergy in MeV
  */

int scatterRspace::getSmatrix(double Ecm0, double Elab0)
{
  cout << "ECM: " << Ecm0 << " Elab0: " << Elab0 << endl;
  Pot->setEnergy(Ecm0);
  newEnergy(Ecm0,Elab0);


  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave;

  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];


  int l = 0;
  lStop = 20;//originally: -1
    //for (;;) // loop over orbital angular momentum
    for (int ii = 0; ii<21 ; ++ii) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          double j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;

	      Pot->setAM(l,j); // set angular momentum of potential

           
          //load in values for wavefunction array and determine where 
          //to start. for very high orbital angular momentum, 
	  //its best not to start to close to the center as the forth 
          //order differential equation solver cannot not handle the 
	  // high derivatives of the centrifugal potential.
          int i = 0;
	    for (;;)
	    {
               y[i] = complex<double>(pow(x[i],l+1),0.);
             // we changed the starting point for integrating 
             // schrod. eq. so we are able to go higher in L"(for MSU work
             // but it seems that does not effect the results
            //   y[i] = complex<double>(pow(x[i],l+1)*1.e-100,0.);
               if (i > 0 &&  x[i] > (double)(l*(l+1))*.00163)break;
               i++;
      //         cout<<"x = "<<x[i]<<"y = "<<y[i]<<endl;
	    } 
	    int istart  = i;
         //   cout<<"istart is = "<<istart<<endl;     

	    integrateWaveFunction(istart,l);

	  complex<double> y_end = y[Nnumerov];
          complex<double> dydr_end = (y[Nnumerov+1] - y[Nnumerov-1])/2./dx; 

//       cout<<"dx = "<< dx<<"yend = "<<y_end<<"\t"<<"dydr = "<<dydr_end<< endl;
	  //rescale in case the numbers are too big, 
	  //this happend for large l with large beta sometimes.
	  //rescaling stops overflow latter in the code
	  complex<double> gg = y_end;
          complex<double> ff = dydr_end;
          dydr_end /= y_end.real();
          y_end /= y_end.real();

    //   cout<<"dx = "<< dx<<"yend = "<<y_end<<"\t"<<"dydr = "<<dydr_end<< endl;



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
	  double SReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double SImag = 2.*(AA*BB+CC*DD)/denominator;

          int up = (plusMinus+1)/2;  //=0 form down and =1 for up
	  phaseShift[l][up] = atan2(SImag,SReal);
          if (phaseShift[l][up] < 0.) phaseShift[l][up] += 2.*pi;
          phaseShift[l][up] /=2.;
         
	  S2[l][up] = pow(SReal,2) + pow(SImag,2);
          S[l][up] = complex<double>(SReal,SImag);
          cout << l << " " << j << " " << Ecm << " " << SReal << " " << SImag << " " << phaseShift[l][up] << endl;

       //	cout<<Ecm<<" "<< "Elab "<<Elab << "  "<<l <<" "<<j<<" "<<SReal<<"  "<<SImag <<endl;

   	//  if (Ecm >50. )
   	  //   {
          //  	cout<<Ecm<<" "<<  l <<" "<<j<<" "<<SReal<<"  "<<SImag <<endl;
          //  } 
   	//  if (l==0 & j==.5 )
   	    // {
   	  // if (Ecm>19  ) 
//               cout<<"Elab = " << Elab<< " Energ = " <<Ecm<<" L = " << l 
  //                 <<" J =  "<<j<<" S real =  "<<SReal<<" S imag =  "<<SImag
    //              <<" S^2  "<<S2[l][up]<<" "<<" phase shift up = "<< phaseShift[l][up]<<endl;
          //  } 
//cout <<y_end.real()<<"\t"<<y_end.imag()<<"\t"<<dydr_end.imag()<<"\t"<<dydr_end.real()<<endl;
	}
      //check to see if we have included enough l-waves
      if (1.-S2[l][1] < 0.001 && 1.-S2[l][0] < 0.001 )
  //    if (1.-S2[l][1] < 0.0000000001 && 1.-S2[l][0] < 0.0000000001 )
	{
          lStop = l;
	  break;
	}
      if (isnan(S2[l][1]) || isnan(S2[l][0]))
	 {
	   lStop = l-1;
           cout<<"NANANANA"<<endl;  
	   break;
	 } 
     // cout<<Ecm0<<"\t"<<l<<endl;
      l++;
      if (l > lMax) break;
    }

      if (lStop == -1) cout << "increase lMax" << endl;
      return 1;
}

//**********************************************************
  /**
   * integrates the nonlocal Schrodinger equation.
   * The wave function is found in the arrays x and y
   \param istart starting position in the x and y arrays for integration
   \param l - orbital angular momentum quantum number
  */

void scatterRspace::integrateWaveFunction(int istart, int l)

  {
    nonLocalArray0(l);
    solveNumerov(istart);


    //if potential is local, then we are finished
    if (Pot->beta_max == 0.)return;



    //correct the local wavefunction for nonlocality
    //the correction at r=0, must be unity, so I scaled the usual correction
    //to make this so.
    Pot->HartreeFock.setEquivLocalDepth(1.- sqrt(mstar));
    for (int i=0;i<=Nnumerov+1;i++) y[i] *=
      (1. +  Pot->HartreeFock.potentialEquivLocal(x[i]))/sqrt(mstar); 

  //  for (int itry = 1;itry<20;itry++)
    for (int itry = 1;itry<5;itry++)
      {
	nonLocalArray1(l);
        solveNumerov(istart);

      }

  }
//*********************************************************
  /**
   * sets the source and F arrays for the the first
   * iteration of the wavefunction for the nonlocal calculation
   \param l orbital angular momentum quantum number
   */
void scatterRspace::nonLocalArray0(int l)
{
  
  // estimate equivalent local volume potential
  // with contributions from hartree fock and the voume and it 
  // dispersive correction
  coulDis = 1.73*Z*Zp/Pot->Rc;

  
 // complex<double> U (-Pot->HartreeFock.Vvol-Pot->VolumeAbove.Vvol,
   //                                        -Pot->VolumeAbove.Wvol);
  complex<double> U (-Pot->HartreeFock.Vvol-Pot->VolumeAbove.Vvol,
                                           -Pot->VolumeAbove.Wvol);

  complex<double> k2;
  for (int i=0;i<1000;i++)
    {
      k2 = (complex<double>(Ecm-coulDis) - U)*kconstant*mu;

      complex <double> Unew = 
        - Pot->HartreeFock.Vvol/(1.+Pot->HartreeFock.A) 
        * exp(-k2*pow(Pot->HartreeFock.beta0,2)/4.) 
        - Pot->HartreeFock.A*Pot->HartreeFock.Vvol/(1.+Pot->HartreeFock.A)
        * exp(-k2*pow(Pot->HartreeFock.beta1,2)/4.) 
        - complex<double>(Pot->VolumeAbove.Vvol,Pot->VolumeAbove.Wvol)
        * exp(-k2*pow(Pot->VolumeAbove.beta,2)/4.);
        - complex<double>(Pot->VolumeBelow.Vvol,Pot->VolumeBelow.Wvol)
          *exp(-k2*pow(Pot->VolumeBelow.beta,2)/4.);

      U = ( U + Unew ) / 2.;

    }

  Pot->HartreeFock.setEquivLocalDepth(-U.real());
  Pot->VolumeAbove.setEquivLocalDepth(-U.imag());

  //effective mass relative to nucleon mass
  mstar = real(
    20.9205*k2/(complex<double>( Ecm + Pot->HartreeFock.Vvol + 
				Pot->VolumeAbove.Vvol + Pot-> VolumeBelow.Vvol - coulDis,
                Pot->VolumeAbove.Wvol ) ) );

    double factorI1 = exp(-pow(Pot->VolumeAbove.beta*Kwave/2.,2));
    double factorI2 = exp(-pow(Pot->VolumeBelow.beta*Kwave/2.,2));

    for (int i=0;i<N+2;i++)
      {
       double fr = Pot->HartreeFock.potentialEquivLocal(x[i]);

       //double two = Pot->Volume.DispersiveCorrection(x[i])*factorI;
       double two = 0.;//Pot->Volume.DispersiveCorrection(x[i])*factorI;
                       //now incorporated in effect local ponetial of 
                       //Hartree Fock
       double three = Pot->SurfaceAbove.DispersiveCorrection(x[i])*factorI1;
       double four = Pot->SurfaceBelow.DispersiveCorrection(x[i])*factorI2;
       fr += two + three + four;


       fr += Pot->coulomb(x[i]);
       fr += Pot->HartreeFock.potentialSurface(x[i]);
       // next the spin-orbit potential
       fr += Pot->SpinOrbit.RealPotential(x[i]);


       fr *= gammaRel*muhbar;

       if (l > 0) fr += (double)(l*(l+1))/pow(x[i],2);

       //if (Ecm > 190 && l == 19 ) fr = (double)(l*(l+1))/pow(x[i],2);
       fr = -(Kwave2 - fr);

       
       double fi = Pot->VolumeAbove.potentialEquivLocal(x[i])
                 + Pot->SurfaceAbove.ImaginaryPot(x[i]) * factorI1;
       /*	  
       double fi = (Pot->Volume.ImaginaryPot(x[i])+
                    Pot->Surface.ImaginaryPot(x[i]))
                   *factorI;
       */
       // next the spin-orbit potential
       fi += Pot->SpinOrbit.ImaginaryPotential(x[i]);
       //if (Ecm > 190 && l == 19) fi = 0.;

       fi *= gammaRel*muhbar;

       f[i] = complex<double>(fr,fi);

       source[i] = complex<double>(0.,0.);

      }

}
//*****************************************************************
  /**
   * Sets the source and F arrays for subsequent iterations
   * of the wavefunction in the nonlocal calculation
   \param l orbital angular monentum quantum number
  */
void scatterRspace::nonLocalArray1(int l)
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

         //if (i < Nnumerov+1) slope = (y[i+1]- y[i-1])/2./dx;
         //else slope = (3.*y[i] + y[i-2] - 4.*y[i-1])/2./dx;

	 if (i < Nnumerov + 2) slope = (-y[i-1]/3. 
                                 -y[i]/2. + y[i+1] - y[i+2]/6.)/dx;
         else if (i < Nnumerov + 1) slope = (y[i+1]- y[i-1])/2./dx;
         else slope = (3.*y[i] + y[i-2] - 4.*y[i-1])/2./dx;

         complex<double> stuff= y[i]/slope;
      
         F = exp(-100.*pow(abs(stuff),2));
         stuff = pow(x[i],l+1)/y[i] -1.;

         F *= 1. - exp(-100.*pow(abs(stuff),2));
	}

      // now do the convolution
      complex<double> v(0.,0.);
      double r1 = x[i];
      for (int j=0;j<=Nnumerov+1;j++)
	{
	  double r2 = x[j];
          if (r1 - r2 > 3.*Pot->beta_max) continue;
          if (r2 - r1 > 3.*Pot->beta_max) break;
	  //double ddr = abs(r1-r2)/3.;


        // ddr is calculated in 'nonlocalPart'
        v += dx * Pot->nonlocalPart( r1, r2 ) * y[j];

/*
	  if (Pot->HartreeFock.beta0 > 0. && ddr < Pot->HartreeFock.beta0)
	    v += dx*Pot->angleIntegration(r1,r2,Pot->HartreeFock.beta0,l)*y[j]*
	      complex<double>(Pot->HartreeFock.U0(r1,r2),0.);


	  if (Pot->HartreeFock.beta1 > 0. && ddr < Pot->HartreeFock.beta1 
                 && Pot->HartreeFock.A != 0)
	    v += dx*Pot->angleIntegration(r1,r2,Pot->HartreeFock.beta1,l)*y[j]*
	      complex<double>(Pot->HartreeFock.U1(r1,r2),0.);


	  if (Pot->Volume.beta > 0. && ddr < Pot->Volume.beta)
	    {
	      complex<double> fs = Pot->Surface.U(r1,r2);
	      complex<double> fv = Pot->Volume.U(r1,r2);

	      v += dx*Pot->angleIntegration(r1,r2,Pot->Volume.beta,l)*y[j]*(fv+fs);
	    }
*/

	}


      source[i] = muhbar*F*v*gammaRel;
      complex <double> vloc = muhbar*(1.-F)/y[i]*v*gammaRel;

      // now add in the remainding local parts of the potential
      complex<double> VlocalPart = Pot->localPart(x[i]);
      VlocalPart *= gammaRel*muhbar;
      VlocalPart += complex<double>((double)(l*(l+1))/pow(x[i],2) - Kwave2);

      f[i] = vloc + VlocalPart;



    }


}
//****************************************************************
  /**
   * returns the elastic center-of-mass differntial cross section.
   * the function getSmatrix() must be executed first
   \param theta center-of-mass angle in radians 
  */
double scatterRspace::DifferentialXsection(double theta)
{
  complex<double> A(0.,0.);
  complex<double> B(0.,0.);
  complex<double>tempA;
  complex<double>tempB;

  int ll = lMax;
  legendre Poly(ll);
   
   
  // start with l=0
  // eta is the S matrix

  A = S[0][1];
  A -= 1.;

  if (Zp == 1.) A *=  exp(complex<double>(0.,2.*Sigma[0])); 
  A *= Poly.LegendreP0(0,theta);


  //now the other l's

  for (int i=1;i<=lStop;i++)
    {

      //spin nonflip
      int l = i;

      tempA = (double)(l+1)*S[l][1] + (double)l*S[l][0] ;


      tempA -=  complex<double>(2*l+1);
      tempA *= Poly.LegendreP0(l,theta);
      

      //spin flip
      tempB = (S[l][1]-S[l][0])*Poly.LegendreP1(l,theta);

      if (Zp == 1.) // include Coulomb phase shift
	{
	  tempA *= exp(complex<double>(0.,2.*Sigma[l]));
          tempB *= exp(complex<double>(0.,2.*Sigma[l]));

	}

      A += tempA;
      B += tempB;
    }



  A /= complex<double>(0.,2.*Kwave);
  B /= complex<double>(0.,2.*Kwave);  



  if (Zp == 1)
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
   * the function getSmatrix() must be executed first
   */
double scatterRspace::AbsorptionXsection()
{
  //l==0 contribution
  SigmaAb[0] = (1. - S2[0][1]) *konst;
  SigmaAb_pos[0] = SigmaAb[0];
  SigmaAb_neg[0] = 0.; 
  double tot = SigmaAb[0];
  //contribution from higher l waves
 


  for (int i=1;i<=lStop;i++)
    {
      int l = i;
      int parity = 1-2*(i%2);
      double up = (double)(l+1)*(1.-S2[l][1])*konst;
      double down = (double)l*(1.-S2[l][0])*konst;
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
   * returns the total cross section, infinite for protons (don't try),
   * finite for neutrons in mb
   * The function getSmatrix() must be called first
   */
double scatterRspace::TotXsection()
{
  //l==0 contribution
  double tot = 1.-real(S[0][1]);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*(1.-real(S[i][1])) +
      (double)i*(1.-real(S[i][0]));
    }
  return tot*konst*2.;
  
} 
//*************************************************************
  /**
   * return the total Elastic cross section, 
   * infinite for protons (don't try), finite for neutrons in mb
   * The function getSmatrix() must be called first
   */

double scatterRspace::ElasticXsection()
{
  complex<double> one(1.,0.);
  //l==0 contribution
  double tot = pow(abs(one-S[0][1]),2);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*pow(abs(one-S[i][1]),2) + 
	(double)i*pow(abs(one-S[i][0]),2);
    }
  return tot*konst;
  
} 
//**************************************************************
  /**
   * returns the Rutherford differential cross section in mb
\param theta is the center-of-mass scattering angle in radians
  */
double scatterRspace::Rutherford(double theta)
{
  return 10.* pow(gamma,2)/4./Kwave2/pow(sin(theta/2.),4);
}
//**************************************************************
  /**
   *returns the transmission coefficients
   *the function integrateWave() must be executed first
   \param l is the orbital angular momentum of the nucleon
   \param j is the total angular momentum of the nucleon
   */
double scatterRspace::TransCoef(int l, double j)
{
  if (l > lStop) return 0.;
  int jj = (int)j;
  if (jj == l)return (1.-S2[l][1]);
  else return 1. - S2[l][0];
}
