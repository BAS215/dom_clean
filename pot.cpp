#include "pot.h"

double const pot::pi = acos(-1.);

/**
 * initialized the calss for a new target projectile
\param Z0 - target proton number
\param Zp0 - projectile proton number
*/
void pot::init( double Z0, double Zp0, double A0, int readCoulomb0 )
{
  Z = Z0;
  Zp = Zp0;
  A = A0;
  readCoulomb =  readCoulomb0;
 // if ( readCoulomb == 1 ) read_chd_from_file( A, Z );
  if ( readCoulomb == 1 || 2 ) read_chd_from_file( A, Z );
//  read_chd_from_file( A, Z );
}
//******************************************
  /**
   * Constructor
   \param typeN0 0= sqrt(f(r1)*f(r2))   1=(r1+r2)/2
   */
pot::pot(int typeN0)
{
  typeN = typeN0;
  setTypeN(typeN);
  Ecm = 9999; // just to make sure that it is initialized.
}
//*********************************************************
  /**
   *specifies the type of nonlocality
   \param typeN0 0= sqrt(f(r1)*f(r2))   1=(r1+r2)/2
   */
void pot::setTypeN(int typeN0)
{
  typeN = typeN0;
  HartreeFock.setTypeN(typeN);
  SurfaceAbove.setTypeN(typeN);
  SurfaceBelow.setTypeN(typeN);
  VolumeAbove.setTypeN(typeN);
  VolumeBelow.setTypeN(typeN);
}
//***************************************************
  /**
   * sets the energy and calculates the strengths of the real and imaginary 
   * components of the potential
   \param Ecm0 - center-of-mass energy in MeV
  */
void pot::setEnergy(double Ecm0)
{
  Ecm = Ecm0;
  SurfaceAbove.setEnergy(Ecm);
  SurfaceBelow.setEnergy(Ecm);
  VolumeAbove.setEnergy(Ecm);
  VolumeBelow.setEnergy(Ecm);
  SpinOrbit.setEnergy(Ecm);
}
//******************************************************
  /**
   * returns the local potential.
   * setEnergy(double) and setAM(int,double) must be run before hand 
   * to specify the energy and angular momenta.
   \param r is radius in fm
  */
complex<double>pot::potential(double r)
{
// double vreal =  coulomb(r) + coulomb_exp_screen(r) + HartreeFock.potential(r);
  double vreal = coulomb(r) + HartreeFock.potential(r);
  complex<double> out = complex<double>(vreal,0.);
  out += SpinOrbit.potential( r ) + SurfaceAbove.potential( r ) 
       + SurfaceBelow.potential( r ) + VolumeAbove.potential( r )
       + VolumeBelow.potential( r );
  return out;
}
//*******************************************************
  /**
   * returns the nonlocal potential
   * setEnergy(double) and setAM(int,double) must be run beforehand 
   * to specify the enrgy and angular momentum 
   \param r1 is the magnitude of the first radius
  \param r2 is the magnitude of the second radius
  \param dr is the magnitude of the separation between these radii
  */
/*
complex<double>pot::potential(double r1, double r2, double dr)
{
  complex<double> out;
  //start with the local part
  if (dr == 0) 
    {
      out = complex<double>(coulomb(r1),0.) + SpinOrbit.potential(r1);
    }
  out += complex<double>(HartreeFock.potential(r1,r2,dr),0.);
  out += SurfaceAbove.potential(r1,r2,dr) + SurfaceBelow.potential(r1,r2,dr)
       + VolumeAbove.potential(r1,r2,dr) + VolumeBelow.potential(r1,r2,dr);

  return out;
}*/
//********************************************************
  /**
   * returns the local potential for a given energy and radius
   * setAM(int,double) must be rum beforehand 
   * to specify the angular momenta.
   \param r is the radius in fm
   \param Ecm is the center-of-mass energy in MeV
  */
complex<double>pot::potentialE(double r, double Ecm)
{
  setEnergy(Ecm);
  return potential(r);
}
//*********************************************************
  /**
   * returns the nonlocal potential.
   * setAM(int,double) must be run before hand to specify the 
   * l,j values.
   \param r1 - magnitude of first position vector in fm
   \param r2 - magnitude of second position vector in fm
   \param dr - magnitude of different between two position vectors in fm
   \param Ecm - center-of-mass energy in MeV 
  */
/*
complex<double> pot::potentialE(double r1, double r2, double dr, double Ecm)
{
  setEnergy(Ecm);
  return potential(r1,r2,dr);
}*/
//**********************************************************
  /**
   * sets the l and J values
   \param l0 - orbital AM
   \param j0 - total AM 
  */
void pot::setAM(int l0, double j0)
{
  SpinOrbit.setAM(l0,j0);
  l = l0;
  j = j0;
}
double pot::coulomb_point_screen( double r ) {

    double RNew = 25.;
    double PowerNew = 4.;
    double exp_term = exp(-pow((r/RNew),PowerNew));

    return exp_term*e2*(Z)/r;
}
double pot::coulomb_homogenous( double r ) {

    if (Zp == 0.) return 0.;
    if (r > Rc) return e2*(Z)/r;
    else return  e2*(Z)/2./Rc*(3.-pow(r/Rc,2));
}
double pot::coulomb_exp( double r ) {

    if (rmesh_chd.empty() ) {
        cout<<readCoulomb<<endl;
        cout << "rmesh_chd is empty" << std::endl;
        std::abort();   
    }

    if (Zp == 0.) return 0.;
    double deltarp=rmesh_chd[1]-rmesh_chd[0];	
    double sum=0;
    for (unsigned int i=0; i < rmesh_chd.size() ; ++i) {
	if ( r > rmesh_chd[i]) sum=sum+exp_charge_density[i] * std::pow(rmesh_chd[i] , 2) / r * deltarp; 
	else sum=sum+exp_charge_density[i] * rmesh_chd[i] * deltarp; 

    }
return e2 * 4. * M_PI * sum;
}
//*************************************************************
//********Scrren Coulumb added to make the FB transform possible
//********PhysRevC 71, 054005(2005) A. Deltuva ,..***********
//*************************************************************
double pot::coulomb_exp_screen( double r ) {

    double RNew = 5.0;
    double PowerNew = 4.;
    double exp_term = exp(-pow((r/RNew),PowerNew));
    if (rmesh_chd.empty() ) {
      //  cout << "rmesh_chd is empty" << std::endl;
        cout << "rmesh_chd is empty" << std::endl;
        std::abort();   
    }
 //  cout<<"Whats up screen coulomb"<<endl;
    if (Zp == 0.) return 0.;
    double deltarp=rmesh_chd[1]-rmesh_chd[0];	
    double sum=0;
    for (unsigned int i=0; i < rmesh_chd.size() ; ++i) {
    //    cout<<rmesh_chd[i]<<" "<<  exp_charge_density[i]<<endl; 
	if ( r > rmesh_chd[i]) sum=sum+exp_charge_density[i] * std::pow(rmesh_chd[i] , 2) / r * deltarp; 
	else sum=sum+exp_charge_density[i] * rmesh_chd[i] * deltarp; 

    }
return e2 * 4. * M_PI * sum * exp_term;
}
//*************************************************************
//*************************************************************
//*************************************************************
//**:***********************************************************
//
  /**
   * returns the Coulomb potential in MeV
   \param r = radius in fm
  */
double pot::coulomb(double r)
{

    if ( readCoulomb == 0 ) return coulomb_homogenous( r ); 
   // if ( readCoulomb == 0 ) return coulomb_exp_screen( r ); 
    else if ( readCoulomb == 1 )  return coulomb_exp( r );
    else if ( readCoulomb == 2 ) return coulomb_exp_screen( r );
   // else if ( readCoulomb == 1 ) return coulomb_exp( r );
    else {

        std::cout << "In Function pot::coulomb: " << std::endl;
        std::cout << "Invalid value for readCoulomb." << std::endl;
        std::cout << "readCoulomb = " << readCoulomb << std::endl;
        std::abort();

  }
  
}

//*************************************************************
  /**
   * load the parametrs defining the DOM potential
   \param Rc0 Coulomb radius in fm
   \param VHFvol - depth of HartreeFock volume term in Mev
   \param VHFsur - strength of Hartree Fock surface term in MeV
   \param RHF - radius of Hartree Fock potential in fm
   \param aHF - diffuseness of Hartree Fock potential in fm
   \param beta_nl_R0 - first nonlocality length in fm
   \param AHF - relative contribution from second nonlocality
   \param beta_nl_R1 - second nonlocaility length in fm
   \param Rsur - radius of imaginary surface potential in fm
   \param asur - diffuseness of imaginary surface potential in fm
   \param Asur - Strength of the imaginary surface potential in MeV
   \param Bsur - paramter for imaginary surface
   \param Csur - paramter for imaginary surface
   \param Dsur - parameter for imaginary surface
\param WstartSur - energy from fermi energy where surface starts (MeV)
\param Efermi - Fermi energy in MeV
\param beta_nl_I - nonlocal length for imaginary potentials in fm
\param Rzero - basic radius for volume imaginary in fm
\param deltaR - radius change for volume imaginary term in fm
\param expR - change of radius eith energy parameter 
\param avol - diffuseness for volume potential
\param Avol - strength of imaginary volume potential in MeV
\param Bvol - parameter for imaginary volume
\param Epvol - energy from Fermi where imaginary Volume starts
\param m - exponent for the imaginary volume potential
\param Asy - switch to include energy-saymmetry term
\param alpha - energy asymmetry parameter
\param Ea - energy asymmetry parameter
\param Rso - radius for spin-orbit potential
\param aso - diffuseness for spin orbit potential
\param Vso - strength of SPin Orbit potrential
\param AWso - parametr for imaginary spin-orbit potential
\param BWso - parameter for imaginary spin-orbit potential
   */
/*
void pot::load(double Rc0, double VHFvol, double VHFsur, double RHF, 
               double aHF, double beta_nl_R0, double AHF, double beta_nl_R1,
               double Rsur, double asur, double AsurAbove, double AsurBelow, 
               double Bsur, double Csur, double Dsur, double WstartSur_A, 
               double WstartSur_B, double Efermi, double beta_nl_I0, 
               double beta_nl_I1, double beta_nl_I1_sur, double beta_nl_I0_sur,
	       double Rzero, double deltaR, double expR,
               double avol, double AvolAbove, double AvolBelow, double Bvol,
               double Epvol, int m, int Asy, double alpha, double Ea_above,
               double Ea_below, double Rso, double aso, double Vso,
               double AWso, double BWso) {

  Rc = Rc0;
  HartreeFock.load( VHFvol, VHFsur, RHF, aHF, beta_nl_R0, AHF, beta_nl_R1);

  SurfaceBelow.load( Rsur, asur, AsurBelow, Bsur, Csur, Dsur, 
                     WstartSur_B, Efermi, beta_nl_I0_sur, 0, 0 );

  SurfaceAbove.load( Rsur, asur, AsurAbove, Bsur, Csur, Dsur, 
                     WstartSur_A, Efermi, beta_nl_I1_sur, 1, 0 );

  //VolumeBelow.load( Rzero, deltaR, expR, avol, AvolBelow ,Bvol, Epvol,
  //                  Efermi, m, Asy, alpha, Ea, beta_nl_I0 );

  VolumeBelow.load( Rzero, deltaR, expR, avol, AvolBelow ,Bvol, Epvol,
                    Efermi, m, Asy, alpha, Ea_below , beta_nl_I0 ,0 );

  VolumeAbove.load( Rzero, deltaR, expR, avol, AvolAbove ,Bvol, Epvol,
                    Efermi, m, Asy, alpha, Ea_above , beta_nl_I1 ,1 );

  SpinOrbit.load( Rso, aso, Vso, Efermi, AWso, BWso );

  //find largest nonlocality length
  beta_max = max(beta_nl_R0,beta_nl_I0);
  if (AHF != 0.) beta_max = max(beta_max,beta_nl_R1);
}*/
//*************************************************************
  /**
   * load the parametrs defining the DOM potential
   \param Rc0 Coulomb radius in fm
   \param VHFvol - depth of HartreeFock volume term in Mev
   \param VHFsur - strength of Hartree Fock surface term in MeV
   \param RHF - radius of Hartree Fock potential in fm
   \param aHF - diffuseness of Hartree Fock potential in fm
   \param beta_nl_R0 - first nonlocality length in fm
   \param AHF - relative contribution from second nonlocality
   \param beta_nl_R1 - second nonlocaility length in fm
   \param Rsur - radius of imaginary surface potential in fm
   \param asur - diffuseness of imaginary surface potential in fm
   \param Asur - Strength of the imaginary surface potential in MeV
   \param Bsur - paramter for imaginary surface
   \param Csur - paramter for imaginary surface
   \param Dsur - parameter for imaginary surface
\param WstartSur - energy from fermi energy where surface starts (MeV)
\param Efermi - Fermi energy in MeV
\param beta_nl_I - nonlocal length for imaginary potentials in fm
\param Rzero - basic radius for volume imaginary in fm
\param deltaR - radius change for volume imaginary term in fm
\param expR - change of radius eith energy parameter 
\param avol - diffuseness for volume potential
\param Avol - strength of imaginary volume potential in MeV
\param Bvol - parameter for imaginary volume
\param Epvol - energy from Fermi where imaginary Volume starts
\param m - exponent for the imaginary volume potential
\param Asy - switch to include energy-saymmetry term
\param alpha - energy asymmetry parameter
\param Ea - energy asymmetry parameter
\param Rso - radius for spin-orbit potential
\param aso - diffuseness for spin orbit potential
\param Vso - strength of SPin Orbit potrential
\param AWso - parametr for imaginary spin-orbit potential
\param BWso - parameter for imaginary spin-orbit potential
   */
void pot::load(double Rc0, double VHFvol, double VHFsur,
	       double RHF, double aHF, double RHFs, double aHFs,
	       double beta_nl_R0, double AHF, double beta_nl_R1,
               double RsurAbove, double RsurBelow, double asurAbove,double asurBelow, double AsurAbove, double AsurBelow, 
               double BsurA, double CsurA, double DsurA, double Bsur, double Csur, double Dsur, double WstartSur_A, 
               double WstartSur_B, double Efermi, double beta_nl_I0, 
               double beta_nl_I1,double beta_nl_I0_sur,double beta_nl_I1_sur,
	       double RzeroAbove, double RzeroBelow,
	       double deltaR, double expR,
               double avolAbove, double avolBelow, double AvolAbove, double AvolBelow,
	       double BvolAbove, double BvolBelow, double EpvolAbove, double EpvolBelow,
	       int m, int Asy, double alpha, double Ea_above,
               double Ea_below, double Rso, double aso, double Vso,
               double AWso, double BWso , double V_wine1 , double R_wine1 , double rho_wine1 ) {

  Rc = Rc0;
  HartreeFock.load( VHFvol, VHFsur, RHF, aHF, RHFs, aHFs, beta_nl_R0, AHF, beta_nl_R1 , V_wine1 , R_wine1, rho_wine1 );

  SurfaceBelow.load( RsurBelow, asurBelow, AsurBelow, Bsur, Csur, Dsur, 
                     WstartSur_B, Efermi, beta_nl_I0_sur, 0, 0 );

  SurfaceAbove.load( RsurAbove, asurAbove, AsurAbove, BsurA, CsurA, DsurA, 
                     WstartSur_A, Efermi, beta_nl_I1_sur, 1, 0 );

  //VolumeBelow.load( Rzero, deltaR, expR, avol, AvolBelow ,Bvol, Epvol,
  //                  Efermi, m, Asy, alpha, Ea, beta_nl_I0 );

  VolumeBelow.load( RzeroBelow, deltaR, expR, avolBelow, AvolBelow ,BvolBelow, EpvolBelow,
                    Efermi, m, Asy, alpha, Ea_below , beta_nl_I0 ,0 );

  VolumeAbove.load( RzeroAbove, deltaR, expR, avolAbove, AvolAbove ,BvolAbove, EpvolAbove,
                    Efermi, m, Asy, alpha, Ea_above , beta_nl_I1 ,1 );

  SpinOrbit.load( Rso, aso, Vso, Efermi, AWso, BWso );

  //find largest nonlocality length
  beta_max = max(beta_nl_R0,beta_nl_I0);
  beta_max = max(beta_max,beta_nl_I0_sur);
  if (AHF != 0.) beta_max = max(beta_max,beta_nl_R1);
}
//*******************************************************
  /**
   * returns the local part of the potential.
   * setEnergy(double) and setAM(int,l) must be run beforehand 
   * to specify the energy and angular momentum
   \param r is the radius in fm
  */


complex<double> pot::localPart(double r)
{
// double vreal =  coulomb_exp_screen(r) ;
 double vreal = coulomb(r) ;
// cout<<r<<" "<<coulomb(r)<<endl;
  if (HartreeFock.beta0 == 0)
    vreal += HartreeFock.potential0(r);
  if (HartreeFock.beta1 == 0 && HartreeFock.A != 0. )
    vreal += HartreeFock.potential1(r);

  if (HartreeFock.R_wine == 0)
    vreal += HartreeFock.BWU_local(r);

  vreal += SpinOrbit.RealPotential(r);
  vreal += HartreeFock.potentialSurface(r);

//COMMENTED out for perry buck calculation
//************************************
  double vimag=  SpinOrbit.ImaginaryPotential(r);
  if (SurfaceAbove.betas == 0)
    {
      vreal += SurfaceAbove.DispersiveCorrection(r);
      vimag += SurfaceAbove.ImaginaryPot(r);
    }

  if (SurfaceBelow.betas == 0)
    {
      vreal += SurfaceBelow.DispersiveCorrection(r);
      vimag += SurfaceBelow.ImaginaryPot(r);
    }

  if (VolumeAbove.beta == 0)
    {
      vreal += VolumeAbove.DispersiveCorrection(r);
      vimag += VolumeAbove.ImaginaryPot(r);
    }

  if (VolumeBelow.beta == 0)
    {
      vreal += VolumeBelow.DispersiveCorrection(r);
      vimag += VolumeBelow.ImaginaryPot(r);
    }
/* 
  double W_PB = 10.187;
  double R_surf_PB = 1.250 * pow(40. , 1./3.); 
  double a_surf_PB = .423; 
  double fact_PB = exp((r-R_surf_PB)/a_surf_PB);
  double  vimag_perrybuck = -4.* fact_PB * W_PB/pow(1+fact_PB,2);
   double vreal_perrybuck = -44.224* 1./(1+exp( (r-1.298*pow(40.0,1./3.))/.614) ); 
   vimag += vimag_perrybuck;
   vreal += vreal_perrybuck;
*/
 // vimag = 0;
 // vreal = 0;
 // vreal = vreal+ coulomb(r);
  return complex<double>(vreal,vimag);
}
//*******************************************************
  /**
   * returns the nonlocal part of the potential.
   * setEnergy(double) and setAM(int,l) must be run beforehand 
   * to specify the energy and angular momentum
   \param r is the radius in fm
  */
complex<double> pot::nonlocalPart( double r1, double r2 )
{
        
	double ddr = abs(r1-r2)/3.;
	//double ddr = -1;

    double v_hf = 0.0; // Hartree-Fock
    double v_Nunes = 0.0; // Hartree-Fock
    complex< double > v_dy( 0.0, 0.0 ); // Dynamic


    if ( HartreeFock.R_wine > 0. && ddr < HartreeFock.R_wine ) { 

     v_hf += angleIntegration( r1, r2, HartreeFock.R_wine, l )*HartreeFock.BWU(r1,r2); 
    }

    if ( HartreeFock.beta0 > 0. && ddr < HartreeFock.beta0 ) { 

//      v_hf += HartreeFock.U_Angular_Integration( r1, r2, HartreeFock.beta0, l);
     v_hf += angleIntegration( r1, r2, HartreeFock.beta0, l )*HartreeFock.U0(r1,r2); 
     v_Nunes += angleIntegration( r1, r2, HartreeFock.beta0, l )*HartreeFock.U0(r1,r2); 
    }

    if ( HartreeFock.beta1 > 0. && ddr < HartreeFock.beta1 && 
         ( HartreeFock.A != 0 ) ) {

//      v_hf += HartreeFock.U_Angular_Integration( r1, r2,HartreeFock.beta1, l);
      v_hf += angleIntegration( r1, r2, HartreeFock.beta1, l )*HartreeFock.U1(r1,r2); 
      v_Nunes += angleIntegration( r1, r2, HartreeFock.beta1, l )*HartreeFock.U1(r1,r2); 
    }


    double gaussAbove = 0;
    if ( VolumeAbove.beta > 0 && ddr < VolumeAbove.beta ) {

    gaussAbove = angleIntegration( r1, r2, VolumeAbove.beta, l ); 
 //      gaussAbove = FullAngleIntegration( r1,  r2, VolumeAbove.beta , l , 1.0  ,  VolumeAbove.a , VolumeAbove.Rzero); 
    }

    double gaussBelow = 0;
    if ( VolumeBelow.beta > 0 && ddr < VolumeBelow.beta ) {

       gaussBelow = angleIntegration( r1, r2, VolumeBelow.beta, l ); 
//      gaussBelow = FullAngleIntegration( r1,  r2, VolumeBelow.beta , l , 1.0  ,  VolumeBelow.a , VolumeBelow.Rzero); 

    }
  
    double gaussAbove_sur = 0;
    if ( SurfaceAbove.betas > 0 && ddr < SurfaceAbove.betas ) {

   gaussAbove_sur = angleIntegration( r1, r2, SurfaceAbove.betas, l ); 
//      gaussAbove_sur = FullAngleIntegration_sur( r1,  r2, SurfaceAbove.betas , l , 1.0  ,  SurfaceAbove.a , SurfaceAbove.R); 
    }

    double gaussBelow_sur = 0;
    if ( SurfaceBelow.betas > 0 && ddr < SurfaceBelow.betas ) {

       gaussBelow_sur = angleIntegration( r1, r2, SurfaceBelow.betas, l ); 
//       gaussBelow_sur = FullAngleIntegration_sur( r1,  r2,  SurfaceBelow.betas , l , 1.0  , SurfaceBelow.a , SurfaceBelow.R);
    }


//    v_dy += SurfaceAbove.Un() * gaussAbove_sur; 
    v_dy += SurfaceAbove.U( r1, r2 ) * gaussAbove_sur; 

//    v_dy += SurfaceBelow.Un() * gaussBelow_sur; 
    v_dy += SurfaceBelow.U( r1, r2 ) * gaussBelow_sur; 

//    v_dy += VolumeAbove.Un() * gaussAbove;
   v_dy += VolumeAbove.U( r1, r2 ) * gaussAbove;

//    v_dy += VolumeBelow.Un() * gaussBelow;
    v_dy += VolumeBelow.U( r1, r2 ) * gaussBelow;

  return v_hf + v_dy;
//  return 0;
 // v_dy = 0.;
 //return v_Nunes + v_dy;
}

double pot::nonlocalIM( double r1, double r2 ) {

	double ddr = abs(r1-r2)/3.;
	//double ddr = -1;

    complex< double > v_dy( 0.0, 0.0 ); // Dynamic

    double gaussAbove_sur = 0;
    if ( SurfaceAbove.betas > 0 && ddr < SurfaceAbove.betas ) {

       gaussAbove_sur = angleIntegration( r1, r2, SurfaceAbove.betas, l ); 
//      gaussAbove_sur = FullAngleIntegration_sur( r1,  r2, SurfaceAbove.betas , l , 1.0  ,  SurfaceAbove.a , SurfaceAbove.R); 
    }

    double gaussBelow_sur = 0;
    if ( SurfaceBelow.betas > 0 && ddr < SurfaceBelow.betas ) {

       gaussBelow_sur = angleIntegration( r1, r2, SurfaceBelow.betas, l ); 
//       gaussBelow_sur = FullAngleIntegration_sur( r1,  r2,  SurfaceBelow.betas , l , 1.0  , SurfaceBelow.a , SurfaceBelow.R);
    }

    double gaussAbove = 0;
    if ( VolumeAbove.beta > 0 && ddr < VolumeAbove.beta ) {

       gaussAbove = angleIntegration( r1, r2, VolumeAbove.beta, l ); 
 //     gaussAbove = FullAngleIntegration( r1,  r2, VolumeAbove.beta , l , 1.0  ,  VolumeAbove.a , VolumeAbove.Rzero); 
    }

    double gaussBelow = 0;
    if ( VolumeBelow.beta > 0 && ddr < VolumeBelow.beta ) {

       gaussBelow = angleIntegration( r1, r2, VolumeBelow.beta, l ); 
//      gaussBelow = FullAngleIntegration( r1,  r2, VolumeBelow.beta , l , 1.0  ,  VolumeBelow.a , VolumeBelow.Rzero); 
    }

    // For now, it is assumed that the volume and surface nonlocalities 
    // are the same
    v_dy += SurfaceAbove.U( r1, r2 ) * gaussAbove_sur; 
    v_dy += SurfaceBelow.U( r1, r2 ) * gaussBelow_sur; 

    v_dy += VolumeAbove.U( r1, r2 ) * gaussAbove;
    v_dy += VolumeBelow.U( r1, r2 ) * gaussBelow;

//    v_dy += SurfaceAbove.Un() * gaussAbove_sur; 

//   v_dy += SurfaceBelow.Un() * gaussBelow_sur; 

  //  v_dy += VolumeAbove.Un() * gaussAbove;

 //   v_dy += VolumeBelow.Un() * gaussBelow;


    return imag( v_dy );
}

double pot::nonlocal_surface_IM( double r1, double r2 ) {

	double ddr = abs(r1-r2)/3.;
	//double ddr = -1;

    complex< double > v_dy( 0.0, 0.0 ); // Dynamic

    double gaussAbove_sur = 0;
    if (  SurfaceAbove.betas > 0 && ddr <  SurfaceAbove.betas ) {

       gaussAbove_sur = angleIntegration( r1, r2, SurfaceAbove.betas, l ); 
//      gaussAbove_sur = FullAngleIntegration_sur( r1,  r2, SurfaceAbove.betas , l , 1.0  ,  SurfaceAbove.a , SurfaceAbove.R); 
    }

    double gaussBelow_sur = 0;
    if ( SurfaceBelow.betas > 0 && ddr < SurfaceBelow.betas ) {

       gaussBelow_sur = angleIntegration( r1, r2, SurfaceBelow.betas, l ); 
//       gaussBelow_sur = FullAngleIntegration_sur( r1,  r2,  SurfaceBelow.betas , l , 1.0  , SurfaceBelow.a , SurfaceBelow.R);
    }


    // For now, it is assumed that the volume and surface nonlocalities 
    // are NOT the same
    v_dy += SurfaceAbove.U( r1, r2 ) * gaussAbove_sur; 
     v_dy += SurfaceBelow.U( r1, r2 ) * gaussBelow_sur; 


 //   v_dy += SurfaceAbove.Un() * gaussAbove_sur; 

//    v_dy += SurfaceBelow.Un() * gaussBelow_sur; 

    return imag( v_dy );

}

double pot::nonlocalHF( double r1, double r2 )
{
    
	double ddr = abs(r1-r2)/3.;
	//double ddr = -1;

    double v_hf = 0.0; // Hartree-Fock

    if ( HartreeFock.beta0 > 0. && ddr < HartreeFock.beta0 ) { 

 //       v_hf += HartreeFock.U_Angular_Integration( r1, r2, HartreeFock.beta0, l);
        v_hf +=  HartreeFock.U0( r1, r2 ) * angleIntegration( r1, r2, HartreeFock.beta0, l ); 
    }

  if ( HartreeFock.beta1 > 0. && ddr < HartreeFock.beta1 && 
         ( HartreeFock.A != 0 ) ) {

//        v_hf += HartreeFock.U_Angular_Integration( r1, r2, HartreeFock.beta1, l);
        v_hf += HartreeFock.U1( r1, r2 ) * angleIntegration( r1, r2, HartreeFock.beta1, l );
                
    }


    if ( HartreeFock.R_wine > 0. && ddr < HartreeFock.R_wine ) {

        v_hf +=  HartreeFock.BWU(r1, r2 ) * angleIntegration( r1, r2, HartreeFock.R_wine, l );
              
    }
    return v_hf;
}

double pot::localHF(double r)
{
   // double vreal = coulomb(r)+coulomb_exp_screen(r) ;
    double vreal = coulomb(r);

    if (HartreeFock.beta0 == 0) vreal += HartreeFock.potential0(r);
    if (HartreeFock.beta1 == 0 && HartreeFock.A != 0. ) 
        vreal += HartreeFock.potential1(r);

    vreal += HartreeFock.potentialSurface(r);

    return vreal;
}

// These are needed when calculating the spectroscopic factors
double pot::der_disp_localPart(double r) {

    double vreal = 0;
    vreal += SpinOrbit.DerDispersiveCorrection(r);

    if (SurfaceAbove.betas == 0) {
      vreal += SurfaceAbove.DerDispersiveCorrection(r);
    }

    if (SurfaceBelow.betas == 0) {
      vreal += SurfaceBelow.DerDispersiveCorrection(r);
    }

    if (VolumeAbove.beta == 0) {

      vreal += VolumeAbove.DerDispersiveCorrection(r);
    }

    if (VolumeBelow.beta == 0) {

      vreal += VolumeBelow.DerDispersiveCorrection(r);
    }

    return vreal;
}

double pot::der_disp_nonlocalPart( double r1, double r2 ) {
    
    double vreal = 0;

    if ( SurfaceAbove.betas > 0. ) {

        vreal += angleIntegration( r1, r2, SurfaceAbove.betas, l ) 
               * SurfaceAbove.dU( r1, r2 );
    }

    if ( SurfaceBelow.betas > 0. ) {

        vreal += angleIntegration( r1, r2, SurfaceBelow.betas, l ) 
               * SurfaceBelow.dU( r1, r2 );
    }

    if (VolumeAbove.beta > 0.) {
        vreal += angleIntegration( r1, r2, VolumeAbove.beta, l )
               * VolumeAbove.dU( r1, r2 );
    }

    if (VolumeBelow.beta > 0.) {
        vreal += angleIntegration( r1, r2, VolumeBelow.beta, l )
               * VolumeBelow.dU( r1, r2 );
    }

    return vreal;
}

//************************************************************
  /** 
   * returns the nonlocal potential integrated over the angular coorrdinates
   \param r1 - is the first radius in fm
   \param r2 is the second radius in fm
   * setEnergy(double) and setAM(int,double) must be run sometime before hand
   * to specify the energy and l,j values 
  */
complex<double> pot::potentialL(double r1, double r2)
{
  complex<double> out;
  //first the local part
  if (r1 == r2) out += localPart(r1);


  if (HartreeFock.beta0 > 0.)
    out += angleIntegration(r1,r2,HartreeFock.beta0,l)*
      complex<double>(HartreeFock.U0(r1,r2));

  if (HartreeFock.beta1 > 0. && HartreeFock.A != 0)
    out += angleIntegration(r1,r2,HartreeFock.beta1,l)*
      complex<double>(HartreeFock.U1(r1,r2));

  if (SurfaceAbove.betas > 0.)  
    out += angleIntegration(r1,r2,SurfaceAbove.betas,l)*SurfaceAbove.U(r1,r2);

  if (SurfaceBelow.betas > 0.)  
    out += angleIntegration(r1,r2,SurfaceBelow.betas,l)*SurfaceBelow.U(r1,r2);

  if (VolumeAbove.beta > 0.)
   out += angleIntegration(r1,r2,VolumeAbove.beta,l)*VolumeAbove.U(r1,r2);

  if (VolumeBelow.beta > 0.)
   out += angleIntegration(r1,r2,VolumeBelow.beta,l)*VolumeBelow.U(r1,r2);

   return out;
}
////////////////////////////////////////////////////////////////////////////////////
////**************************************************************************//////
////////////////////////////////////////////////////////////////////////////////////
//Integral of the imaginary part of nonlocal potential: \integral Imag(nonlocalPart(r,r')dv dv'
 
double pot::New_Im_pot_Integral(double E,int rpts ,vector<double> rmesh,double rdelt,  int l,int kind , int trips)
{
  setEnergy(E);
  double sum2=0;
  
  double beta; 
  double R0;
  double a;
  complex<double> U; 
  if (kind == 0 ){
	  if (trips == 0 ){ R0 = SurfaceBelow.R;    beta = SurfaceBelow.betas; a = SurfaceBelow.a; U = SurfaceBelow.Un();} 
 	 
          else if (trips == 2 ){ R0 = SurfaceAbove.R;    beta = SurfaceAbove.betas;  a = SurfaceAbove.a; U = SurfaceAbove.Un(); }

          for ( int ii=0; ii<rpts;++ii){
        	 double sum1=0;  
	         for (int jj=0 ; jj<rpts;++jj){
       		    sum1+= FullAngleIntegration_sur(rmesh[jj],rmesh[ii], beta , l ,  1.0 , a , R0)*rmesh[jj]*rmesh[jj]*rdelt;
                 }
	         sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 

          }
        return sum2 * std::pow(4.0 * M_PI,2) * imag(U);
  }

  else if (kind == 1 ){

 	  if (trips == 1 ){ R0 = VolumeBelow.Rzero; beta = VolumeBelow.beta;   a = VolumeBelow.a; U = VolumeBelow.Un(); }
 
          else if (trips == 3 ){ R0 = VolumeAbove.Rzero; beta = VolumeAbove.beta;  a = VolumeAbove.a;  U = VolumeAbove.Un(); }

	  for ( int ii=0; ii<rpts;++ii){
        	 double sum1=0;  
	         for (int jj=0 ; jj<rpts;++jj){
        	     sum1+= FullAngleIntegration(rmesh[jj],rmesh[ii], beta , l ,  1.0 , a , R0)*rmesh[jj]*rmesh[jj]*rdelt;
         	}
	         sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
          }
         return sum2 * std::pow(4.0 * M_PI,2) * imag(U);
  }

}

////////////////////////////////////////////////////////////////////////////////////
////**************************************************************************//////
////////////////////////////////////////////////////////////////////////////////////
double pot::Old_Im_pot_Integral(double E,int rpts ,vector<double> rmesh,double rdelt,  int l, int trips)
{

  setEnergy(E);
  double sum2=0;
  
  double beta=0.; 

  if (trips == 0 ){  beta = SurfaceBelow.betas; 
	  for ( int ii=0; ii<rpts;++ii){
     		 double sum1=0;  
     		 for (int jj=0 ; jj<rpts;++jj){
          
        		  sum1+= angleIntegration(rmesh[jj],rmesh[ii],beta,l)*
				 real(SurfaceBelow.U(rmesh[ii],rmesh[jj])) * 
				 rmesh[jj] * rmesh[jj] * rdelt;
      		 }
      		sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
  	  }
	 return std::pow(4.0 * M_PI,2) * sum2;}

  else if (trips == 1 ){ beta = VolumeBelow.beta;  
	  for ( int ii=0; ii<rpts;++ii){
                double sum1=0;  
      		for (int jj=0 ; jj<rpts;++jj){
          
          	sum1+= angleIntegration(rmesh[jj],rmesh[ii],beta,l)*
		       real(VolumeBelow.U(rmesh[ii],rmesh[jj])) *
		       rmesh[jj] * rmesh[jj] * rdelt;
      		}
      		sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
 	  }	
 	  return std::pow(4.0 * M_PI,2) * sum2;}

  else if (trips == 2 ){ beta = SurfaceAbove.betas; 
	  for ( int ii=0; ii<rpts;++ii){
     		 double sum1=0;  
	         for (int jj=0 ; jj<rpts;++jj){
          
         	 sum1+= angleIntegration(rmesh[jj],rmesh[ii],beta,l)*
			real(SurfaceAbove.U(rmesh[ii],rmesh[jj])) * 
			rmesh[jj] * rmesh[jj] * rdelt;
      		 }
	         sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
 	  }	
	  return std::pow(4.0 * M_PI,2) * sum2;}

  else if (trips == 3 ){ beta = VolumeAbove.beta; 
	  for ( int ii=0; ii<rpts;++ii){
     		 double sum1=0;  
		 for (int jj=0 ; jj<rpts;++jj){
          
        	 sum1+= angleIntegration(rmesh[jj],rmesh[ii],beta,l)*
			real(VolumeAbove.U(rmesh[ii],rmesh[jj])) *
		        rmesh[jj] * rmesh[jj] * rdelt;
      		 }
	         sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
  	  }
	   return std::pow(4.0 * M_PI,2) * sum2;}

}

////////////////////////////////////////////////////////////////////////////////////
////**************************************************************************//////
////////////////////////////////////////////////////////////////////////////////////
//Integral of Volume imaginary over space J_W
//
//
// double pot::VolumeBelow_Im_integral(double E, double rmax ,int rpts,int l)
// {
//    std::vector<double> rmesh;
//    double rdelt = rmax / rpts;
//    for( int i = 0; i < rpts; ++i ) {
//
//        rmesh.push_back( ( i + 0.5 ) * rdelt );
//       
//    }
//  setEnergy(E);
//  double sum2=0;
//  for ( int ii=0; ii<rpts;++ii){
 //     double sum1=0;  
  //    for (int jj=0 ; jj<rpts;++jj){
 //         
  //        sum1+= angleIntegration(rmesh[jj],rmesh[ii],VolumeBelow.beta,l)*imag(VolumeBelow.U(rmesh[ii],rmesh[jj])) * rmesh[jj] * rmesh[jj] * rdelt;
   //   }
    //  sum2+=sum1 * rmesh[ii] * rmesh[ii] * rdelt; 
 // }
// return sum2;
// }


////////////////////////////////////////////////////////////////////////////////////
////**************************************************************************//////
////////////////////////////////////////////////////////////////////////////////////

//****************************************************************
  /**
   * Calculates a factor obtained after the integration of the 
   * H() part of the nonlocal potential over the angular coordinates
   * returns the integrated Gaussian factor H() multiplied by r1 * r2. 
   \param r1 - first radius  in fm
   \param r2 - second radius in fm
   \param beta - nonlocal length in fm
    \param l - orbital angular momentum quantum number
  */
double pot::FullAngleIntegration(double r1, double r2, double beta , int l, double V ,  double a , double Rzero)
{
 int n_Points = 20;
 vector <double> X(n_Points);
 vector <double> dX(n_Points);
 double r0 = -1.;
 double rn = 1.;
 double Sum = 0.;
 GausLeg(r0 , rn, X, dX);
 for (int i=0 ; i < X.size() ; ++i){
	
       double r_subtract = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * X[i]);

       double r_addition = std::sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * X[i]);

       double H = 1./ ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                  * std::exp( - ( std::pow( ( r_subtract ) / beta, 2 ) ) );

       double Leg_pol = boost::math::legendre_p( l , X[i]);           

       double WS =  -V/(1.+std::exp((r_addition / 2. - Rzero) / a));

       Sum += WS * H * Leg_pol * dX[i];   
 
 }
return 2. * r1 * r2 * Sum;
}

double pot::FullAngleIntegration_sur(double r1, double r2, double beta , int l, double V ,  double a , double Rzero)
{
 int n_Points = 20;
 vector <double> X(n_Points);
 vector <double> dX(n_Points);
 double r0 = -1.;
 double rn = 1.;
 double Sum = 0.;
 GausLeg(r0 , rn, X, dX);
 for (int i=0 ; i < X.size() ; ++i){
	
       double r_subtract = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * X[i]);

       double r_addition = std::sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * X[i]);

       double H = 1./ ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                  * std::exp( - ( std::pow( ( r_subtract ) / beta, 2 ) ) );

       double Leg_pol = boost::math::legendre_p( l , X[i]);           


       double fact = exp((r_addition / 2.-Rzero) / a);
       double d_WS = -4.0 * fact * V / pow((1. + fact),2);

       Sum += d_WS * H * Leg_pol * dX[i];   
 
 }
return 2. * r1 * r2 * Sum;
}

double pot::angleIntegration(double r1, double r2, double beta, int l)
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

/*
// This function uses the modified spherical bessel function of the
// first kind from the Gnu Scientific Library. The function is actually
// scaled by exp(-|x|). Returns the integrated Gaussian factor H
double pot::angleIntegration(double r1, double r2, double beta, int L) {

    double x = 2 * r1 * r2 / std::pow( beta, 2 );

    double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                * gsl_sf_bessel_il_scaled( L, x ); 

    return y;
}
*/

//************************************************************
  /** 
   * returns the nonlocal potential integrated over the angular coorrdinates
   \param r1 - is the first radius in fm
   \param r2 is the second radius in fm
    \param Ecm is the center-of-mass energy in MeV
    * setAM(int,double) must be run sometime beforhand to specify the l,j 
    *values
  */
complex<double> pot::potentialL(double r1, double r2, double Ecm)
{
  setEnergy(Ecm);
  return potentialL(r1,r2);
}

// Wrapper for Bob's potential class.
// This function returns a potential class with 
// all the appropriate initializations.
// >> type indicates form factor used for nonlocal potentials:
//      * 1 is the form of Perey and Buck (average)
//      * 0 is the form of D. Van Neck
// >> mvolume is used in the form for the energy dependence 
//    of the imaginary volume potential. 
// >> AsyVolume is used to specify whether the imaginary volume
// >> potential will have an asymmetry dependence (1 == yes, 0 == no)
// >> tz is the isospin of projectile (-0.5 for neutrons, 0.5 for protons)
// >> Nu holds information about the target (Fermi energy, A, Z, etc.)
// >> p is a struct holding all the parameters
pot get_bobs_pot2( int type, int mvolume, int AsyVolume, double tz, 
		    const NuclearParameters &Nu, const Parameters &p ) {

    double Zp; 
    if ( tz > 0 ) Zp = 1;
    else Zp = 0;
    
    pot Pot( type );
    Pot.init( Nu.Z, Zp, Nu.A, Nu.readCoulomb );
    Pot.load( p.Rc, p.VHFvol, p.VHFsur, p.RHF, p.aHF,  p.RHFs, p.aHFs,p.beta_nl_R0, p.AHF,
              p.beta_nl_R1, p.RsurfaceAbove,p.RsurfaceBelow, p.asurfaceAbove,p.asurfaceBelow, p.AsurfaceAbove, 
              p.AsurfaceBelow, p.BsurfaceA, p.CsurfaceA, p.DsurfaceA, p.Bsurface, p.Csurface, p.Dsurface, 
              Nu.Wgap * p.fGap_A , Nu.Wgap * p.fGap_B , Nu.Ef, p.beta_nl_I0, 
              p.beta_nl_I1, p.beta_nl_I0_sur, p.beta_nl_I1_sur,p.RvolumeAbove,p.RvolumeBelow, p.deltaRvolume, p.expRvolume,
              p.avolumeAbove, p.avolumeBelow, p.AvolumeAbove, p.AvolumeBelow, p.BvolumeAbove, p.BvolumeBelow,
              p.EpvolumeAbove, p.EpvolumeBelow, mvolume, AsyVolume, p.alphaVolume, p.EaVolume_a,
              p.EaVolume_b , p.Rso, p.aso, p.Vso, p.AWso, p.BWso ,  p.V_wine , p.R_wine, p.rho_wine);

    return Pot;
}


// Energy dependent parts of Volume potential. 
complex<double> pot::volumeE( double E ) {

    VolumeAbove.setEnergy( E );
    VolumeBelow.setEnergy( E );

    double vol_Re = VolumeAbove.Vvol + VolumeBelow.Vvol;
    double vol_Im = VolumeAbove.Wvol + VolumeBelow.Wvol;

    return complex< double >( vol_Re, vol_Im );
}

/*
double pot::dispersiveVolumeE_Above( double E ) {

    VolumeAbove.setEnergy( E );

    return VolumeAbove.Vvol; 
}

double pot::dispersiveVolumeE_Below( double E ) {

    VolumeBelow.setEnergy( E );

    return VolumeBelow.Vvol; 
}
*/

double pot::derDispersiveVolumeE( double E ) {

    VolumeAbove.setEnergy( E );
    VolumeBelow.setEnergy( E );

    return VolumeAbove.derVvol + VolumeBelow.derVvol;
}

complex<double> pot::surfaceE( double E ) {

/*
// Gives same result // sjw 05/04/2012
    SurfaceAbove.setEnergy( E );
    SurfaceBelow.setEnergy( E );

    double sur_Re = SurfaceAbove.V + SurfaceBelow.V;
    double sur_Im = SurfaceAbove.W + SurfaceBelow.W;
*/

    double sur_Re = SurfaceAbove.dispersiveE( E ) 
                  + SurfaceBelow.dispersiveE( E );

    double sur_Im = SurfaceAbove.imaginaryE( E )
                  + SurfaceBelow.imaginaryE( E );
    return complex< double >( sur_Re, sur_Im );

}

double pot::derDispersiveSurfaceE( double E ) {

    return SurfaceAbove.derDispersiveE( E ) + SurfaceBelow.derDispersiveE( E );

}

void pot::read_coulomb_from_file( double A0, double Z0 ) {

    if( ( A0 == 20 ) && ( Z0 == 20 ) ) {

        std::list< std::string > input = 
            util::read_commented_file( "ca40_coulomb_potential.dat" ); 

        BOOST_FOREACH( std::string line, input ) {

            std::vector< std::string > vec = util::split( line );

            double r = boost::lexical_cast< double >( vec.at( 0 ) );
            double Ucoul = boost::lexical_cast< double >( vec.at( 1 ) );

            rmesh_coulomb.push_back( r );
            coulomb_from_input.push_back( Ucoul );
        }
        
        double rStart = 0.05;
        double rStop = 12;
        unsigned int rpts = 202;
        if ( ( rmesh_coulomb.front() != rStart ) ||
             ( rmesh_coulomb.back() != rStop ) ||
             ( rmesh_coulomb.size() != rpts ) ) {

            std::cout << "In function 'read_coulomb_from_file':" << std::endl;
            std::cout << "Grid from file is not the same " 
                      << "as the grid used in the code." << std::endl;

            std::cout << "rStart = " << rmesh_coulomb.front() << std::endl;
            std::cout << "rStop = " << rmesh_coulomb.back() << std::endl;
            std::cout << "rpts = " << rmesh_coulomb.size() << std::endl;

            std::abort();
        }

    }
    else {

        std::cout << "No Coulomb Potential Data for A, Z = " 
                  << A0 << ", " << Z0 << std::endl;

    }
}

void pot::read_chd_from_file( double A0, double Z0 ) {

    std::string filename;
    if( ( A0 == 40 ) && ( Z0 == 20 ) ) {

        filename = "exp_chd_ca40_fb.out";
    }
    else if ( ( A0 == 48 ) && ( Z0 == 20 ) ) {

        filename = "exp_chd_ca48_fb.out";
    }
    else if ( ( A0 == 208 ) && ( Z0 == 82 ) ) {

        filename = "exp_chd_pb208_fb.out";
    }

    std::list< std::string > input = util::read_commented_file( filename ); 
    BOOST_FOREACH( std::string line, input ) {

        std::vector< std::string > vec = util::split( line );

        double r = boost::lexical_cast< double >( vec.at( 0 ) );
        double Ucoul = boost::lexical_cast< double >( vec.at( 1 ) );

        rmesh_chd.push_back( r );
        exp_charge_density.push_back( Ucoul );
    }
        
}
