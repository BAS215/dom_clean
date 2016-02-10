#ifndef _pot
#define _pot
#include <complex>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include "surfaceGeneral.h"
#include "volume.h"
#include "hartreeFock.h"
#include "spinOrbit.h"
#include "sphericalB.h"
#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include <new>
#include "io.h"
#include "Gauss_q.h"
/**
 *\brief defines the potential energy
 *
 * class to give all infomation of the potential energy
 * includes the possibility of a nonlocal potential
 */


class pot
{
 public:

  static double const e2 = 1.44;
  static double const kconstant = .048192;
  static double const pi;
  double Z;
  double Zp;
  double A;
  int readCoulomb;
  double Rc;
  double Ecm;
  double beta_max;
  double j; //!< total angular momentum
  int l;//!< orbital angular momentum

  surfaceG SurfaceAbove;
  surfaceG SurfaceBelow;
  volume VolumeAbove;
  volume VolumeBelow;
  hartreeFock HartreeFock;
  spinOrbit SpinOrbit;

  void load(double, double, double, double, double, double, double, double, double, double,
	    double, double, double, double, double, double, double, double, double, double,
 	    double, double, double, double, double, double, double, double, double, double,
	    double, double, double, double, double, double, double, double, double, double, double,
 	    int, int,
	    double, double, double, double, double, double, double, double, double, double, double);
  
 double coulomb_homogenous( double r ); 
  double coulomb_exp( double r ); 
  double coulomb_exp_screen( double r ); 
  double coulomb_point_screen( double r ); 
  double coulomb(double r);
  pot();
  void init( double Z0,  double Zp0, double A0, int readCoulomb0 );

  pot(int typeN0);
  
  complex<double>potential(double r);
  complex<double>potentialE(double r, double Ecm);
  complex<double>potential(double r1, double r2, double dr);
  complex<double>potentialE(double r1, double r2, double dr, double Ecm);
  complex<double>localPart(double r);
  complex<double> nonlocalPart( double r1, double r2 );
  double nonlocalIM( double r1, double r2 ); 
  double nonlocal_surface_IM( double r1, double r2 ); 
  double nonlocalHF( double r1, double r2 );
  double localHF(double r);
  double der_disp_localPart(double r); 
  double der_disp_nonlocalPart( double r1, double r2 );
  complex<double>potentialL(double r1, double r2);
  complex<double>potentialL(double r1, double r2, double Ecm);

  complex<double> volumeE( double E ); 
  double derDispersiveVolumeE( double E ); 
  complex<double> surfaceE( double E ); 
  double derDispersiveSurfaceE( double E ); 

  double angleIntegration(double r1, double r2, double beta, int l);
  double FullAngleIntegration(double , double, double , int ,double , double , double);
  double FullAngleIntegration_sur(double , double, double  , int , double  ,  double  , double );
  double FullAngleIntegrationTest(double r1, double r2, double beta , int l);
  void setEnergy(double Ecm);
  void setTypeN(int typeN0);
  void setAM(int l, double j);
  int typeN;

//  double VolumeBelow_Im_integral(double, double ,int ,int);
//  double VolumeAbove_Im_integral(double, double ,int,int );
 

  double New_Im_pot_Integral(double ,int  ,vector<double> ,double , int ,  int , int );
  double Old_Im_pot_Integral(double ,int  ,vector<double> ,double ,  int , int );
//  double SurfaceBelow_Im_integral(double,int , vector<double>,double, int ,int);
//  double SurfaceAbove_Im_integral(double,int, vector<double>, double, int, int);
 
  private:
  void read_coulomb_from_file( double A0, double Z0 ); 
  void read_chd_from_file( double A0, double Z0 ); 
  std::vector< double > rmesh_coulomb;
  std::vector< double > coulomb_from_input;
  std::vector< double > rmesh_chd;
  std::vector< double > exp_charge_density;



};

pot get_bobs_pot2( int type, int mvolume, int AsyVolume, double tz,
                   const NuclearParameters &Nu, const Parameters &p ); 

#endif
