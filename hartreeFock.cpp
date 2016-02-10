#include "hartreeFock.h"

double const hartreeFock::konst  = 1./pow(acos(-1.),3./2.);
/**
 * constructor
 */
hartreeFock::hartreeFock()
{
  typeN = 0;
}
//*********************************************************
/** constructor
\param typeN0 =0 U function separable, =1 U function like Perey and Buck
*/
hartreeFock::hartreeFock(int typeN0)
{
  typeN = typeN0;
}
//*******************************************************
/** sets the type of U function
\param typeN0 =0 U function separable, =1 U function like Perey and Buck
*/
void hartreeFock::setTypeN(int typeN0)
{
  typeN = typeN0;
}


//***************************************************************************
  /**
   * loads the Hartree Fock call with fixed energy independent potential
   * single gaussian type of nonlocality
    \param Vvol0 is depth of Hartree Fock volume potential in MeV
    \param Vsur0 is depth of Hartree Fock surface potential in MeV
    \param R0 is radius in fm
    \param a0 is diffuseness in fm
    \param beta00 is nonlocality parameter of volume component
  */
void hartreeFock::load(double Vvol0, double Vsur0, double R0v, double a0v, double R0s, double a0s,double beta00)
{
  Rv = R0v;
  Rs = R0s;
  av = a0v;
  as = a0s;
  Vvol = Vvol0;
  Vsur = Vsur0;
  beta0 = beta00;
  A = 0.;
  beta1 = 1.;
 
  volume.init(1.,Rv,av);
  surface.init(Vsur,Rs,as);
}

//***************************************************************************
  /**
   * loads the Hartree Fock call with fixed energy independent potential
   * two gaussian type of nonlocality
    \param Vvol0 is depth of Hartree Fock volume potential in MeV
    \param Vsur0 is depth of Hartree Fock surface potential in MeV
    \param R0 is radius in fm
    \param a0 is diffuseness in fm
    \param beta00 is the first nonlocality parameter of volume component
    \param A0 is the fraction of the second nonlocality parameter
    \param beta10 is the second nonlocality parameter in fm
  */
void hartreeFock::load(double Vvol0, double Vsur0, double R0s, double a0s, double R0v, double a0v,
		       double beta00, double A0, double beta10)
{
  Rv = R0v;
  Rs = R0s;
  av = a0v;
  as = a0s;
  Vvol = Vvol0;
  Vsur = Vsur0;
  beta0 = beta00;
  A = A0;
  beta1 = beta10;
 
  volume.init(1.,Rv,av);
  surface.init(Vsur,Rs,as);
}

//***************************************************************************
  /**
   * loads the Hartree Fock call with fixed energy independent potential
   * single gaussian type of nonlocality
    \param Vvol0 is depth of Hartree Fock volume potential in MeV
    \param Vsur0 is depth of Hartree Fock surface potential in MeV
    \param R0 is radius in fm
    \param a0 is diffuseness in fm
    \param beta00 is nonlocality parameter of volume component
  */
/*
void hartreeFock::load(double Vvol0, double Vsur0, double R0, double a0,
double beta00 , double V_wine1 , double R_wine1, double rho_wine1)
{
  R = R0;
  a = a0;
  Vvol = Vvol0;
  Vsur = Vsur0;
  beta0 = beta00;
  A = 0.;
  beta1 = 1.;
 
  V_wine = V_wine1;
  R_wine = R_wine1;
  rho_wine = rho_wine1;

  volume.init(1.,R,a);
  surface.init(Vsur,R,a);
  bottleW.init(V_wine,R_wine,rho_wine);
}
*/
//***************************************************************************
  /**
   * loads the Hartree Fock call with fixed energy independent potential
   * two gaussian type of nonlocality
    \param Vvol0 is depth of Hartree Fock volume potential in MeV
    \param Vsur0 is depth of Hartree Fock surface potential in MeV
    \param R0 is radius in fm
    \param a0 is diffuseness in fm
    \param beta00 is the first nonlocality parameter of volume component
    \param A0 is the fraction of the second nonlocality parameter
    \param beta10 is the second nonlocality parameter in fm
  */
void hartreeFock::load(double Vvol0, double Vsur0, double R0v, double a0v,double R0s, double a0s,
		       double beta00, double A0, double beta10, double V_wine1 , double R_wine1, double rho_wine1)
{
  Rv = R0v;
  Rs = R0s;
  av = a0v;
  as = a0s;
  Vvol = Vvol0;
  Vsur = Vsur0;
  beta0 = beta00;
  A = A0;
  beta1 = beta10;
 
  V_wine = V_wine1;
  R_wine=R_wine1;
  rho_wine=rho_wine1;

  volume.init(1.,Rv,av);
  surface.init(Vsur,Rs,as);
  bottleW.init(V_wine , 0.0 , rho_wine);
}

//****************************************************************************
  /**
   * returns the local hartree fock potential
   \param r = radius at which potential is calculated
   */
double hartreeFock::potential(double r)
{
  return Vvol*volume.WoodSaxon(r) + surface.DerWoodSaxon(r);
}
//*************************************************************************
  /**
   * returns the nonlocal hartree potential
   \param R1 - first position vector in fm
   \param R2 - second position vector in fm
  */
/* double hartreeFock::potential(double *R1, double *R2)
{
  double dr = 0.;
  double r1 = 0.;
  double r2 = 0.;
  for (int i=0;i<3;i++)
    {
      dr += pow(R1[i]-R2[i],2);
      r1 += pow(R1[i],2);
      r2 += pow(R2[i],2);
    }
  dr = sqrt(dr);
  r1 = sqrt(r1);
  r2 = sqrt(r2);
  return potential(r1,r2,dr);
} 
   \param r1 - magnitude of first position vector in fm
   \param r2 - magnitude of second position vector in fm
    \param dr - magnitude of difference between two position vectors in fm
 

double hartreeFock::potential(double r1, double r2, double dr)
{

  double out;
  if (A == 0.) out = H(dr,beta0);
  else out  = (H(dr,beta0)+A*H(dr,beta1))/(1.+A);

  out *= U(r1,r2);

  if (dr == 0.) out += surface.DerWoodSaxon(r1);
  return out;
  } */
//****************************************************
  /**
   * returns the U function of the Perey and Buck Eq (3)
   * for their nonlocal volume potential
   \param r1 - first radius in fm
   \param r2 - second readius in fm 
   */
double hartreeFock::BWU_local(double r)
{
     return -bottleW.Gauss(r);
}
double hartreeFock::BWU(double r1, double r2)
{
     return -bottleW.Gauss((r1+r2)/2.0);
}
double hartreeFock::U(double r1, double r2)
{
  if (typeN == 0)
     return -Vvol*sqrt((volume.WoodSaxon(r1) - bottleW.Gauss(r1)) * (volume.WoodSaxon(r2) - bottleW.Gauss(r2)));
  else if (typeN == 1)
  // return  Vvol*(volume.WoodSaxon((r1+r2)/2.) - V_wine * std::exp(-std::pow(((r1+r2-R_wine)/(2.0 * rho_wine)),2)));
   // return  Vvol*(volume.WoodSaxon((r1+r2)/2.) - bottleW.Gauss((r1+r2)/2.0));
    return  Vvol*(volume.WoodSaxon((r1+r2)/2.));

  else 
    {
    cout << "typeN = " << typeN << " not known" << endl;
    abort();
    }
}
//*******************************************************

//*****************************************************
  /**
   * if beta0 and beta1 are not equal, then we can consider two
   * volume potentials, this gives the U function value for the first
   \param r1 - first radius in fm
   \param r2 - second readius in fm 
   */
double hartreeFock::U0(double r1, double r2)
{
  return U(r1,r2)/(1.+A);
}
//*****************************************************
  /**
   * if beta0 and beta1 are not equal, then we can consider two
   * volume potentials, this gives the U function value for the second
   \param r1 - first radius in fm
   \param r2 - second readius in fm 
   */
double hartreeFock::U1(double r1, double r2)
{
  return A*U(r1,r2)/(1.+A);
}


//**********************************************************************************************
// Full Integral of HF potential
// as a function of r1, r2 over the angle, the whole U, (I mean U1 and U0) contributions are also taken into account
// (the former  approximation that was : vec_r1+vec_r2 = r1+r2 , )
//this integral over the angel is used directly in pot.cpp, where the Interal over angle is needed
//**********************************************************************************************
double hartreeFock::U_Angular_Integration(double r1, double r2,double beta ,int l)
{
 int n_Points = 20;
 vector <double> X(n_Points);
 vector <double> dX(n_Points);
 double r0 = -1.;
 double rn = 1.;
 double Sum = 0.;
 GausLeg(r0 , rn, X, dX);
 if (beta==beta0){ 
  for (int i=0 ; i < X.size() ; ++i){
	
       double r_subtract = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * X[i]);

       double r_addition = std::sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * X[i]);

       double H0 = 1./ ( std::pow( beta0, 3 ) * std::sqrt( M_PI ) ) 
                  * std::exp( - ( std::pow( ( r_subtract ) / beta0, 2 ) ) );

       double Leg_pol = boost::math::legendre_p( l , X[i]);           

       double WS = (volume.WoodSaxon((r_addition)/2.) - bottleW.Gauss((r_addition)/2.0));

       Sum += WS * H0 * Leg_pol * dX[i];   
 
 }
  return 2. * r1 * r2 * Vvol * Sum/(1+A) ; }
 if (beta==beta1){ 
  for (int i=0 ; i < X.size() ; ++i){
	
       double r_subtract = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * X[i]);

       double r_addition = std::sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * X[i]);

       double H1 = 1./ ( std::pow( beta1, 3 ) * std::sqrt( M_PI ) ) 
                  * std::exp( - ( std::pow( ( r_subtract ) / beta1, 2 ) ) );
       
       double Leg_pol = boost::math::legendre_p( l , X[i]);           

       double WS = (volume.WoodSaxon((r_addition)/2.) - bottleW.Gauss((r_addition)/2.0));

       Sum += WS * H1 * Leg_pol * dX[i];   
 
 }
 return 2. * r1 * r2 * Vvol * Sum * A/(1+A);}
}

//****************************************************
/** 
 * returns the H function of Perey and Buck Eq(3)
 * for their nonlocal potential
 */
double hartreeFock::H(double dr, double beta)
{
  if (beta == 0)
    {
      if (dr == 0.) return 1.;
      else return 0.;
    }
  return konst*exp(-pow(dr/beta,2))/pow(beta,3);
} 
//****************************************************
/**
 * sets the equivalent local potential as per 
 * Perey and Buck
\param VvolE0 is depth of equivalent local volume potential
*/
void hartreeFock::setEquivLocalDepth(double VvolE0)
{
  VvolE = VvolE0;
}
//****************************************************
  /**
   * return the equivalent local potential as
   * per Perey and Buck
    \param r - radius in fm
   */
double hartreeFock::potentialEquivLocal(double r)
{
  return VvolE*volume.WoodSaxon(r);
}
//*******************************************
double hartreeFock::potentialSurface(double r)
{
  return  surface.DerWoodSaxon(r);
  //return  surface.WoodSaxon(r);
}
//**********************************************
double hartreeFock::potential0(double r)
{
  return Vvol*volume.WoodSaxon(r)/(1.+A);
}
//**********************************************
//double hartreeFock::potential1(double r)
//{
//  return A*Vvol*volume.WoodSaxon(r)/(1.+A);
//}
double hartreeFock::potential1(double r)
{
  return Vvol*volume.WoodSaxon(r);
}
