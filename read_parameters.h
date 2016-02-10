#ifndef _READ_PARAMETERS_H_
#define _READ_PARAMETERS_H_


#include <fstream>
#include <iostream>
#include <string>
#include <ios>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "NuclearParameters.h"

struct Raw_Parameters {

  Raw_Parameters() {}

  double rc0;
  double rc1;
  double VHFvol;
  double VHFsur;
  double beta_nl_R0;
  double beta_nl_R1;
  double beta_nl_I0;
  double beta_nl_I1;
  double beta_nl_I0_sur;
  double beta_nl_I1_sur;
  double AHF;
  double rHF0;
  double rHF0s;
  double rHF1;
  double aHF0;
  double aHF0s;
  double aHF1;
  double fGap_A;
  double fGap_B;
  double BsurfaceA;
  double CsurfaceA;
  double DsurfaceA;
  double Bsurface;
  double Csurface;
  double Dsurface;

  double Asurface9P;
  double Asurface9N;
  double Asurface40P_A;
  double Asurface40N_A;
  double Asurface40P_B;
  double Asurface40N_B;
  double Asurface42P;
  double Asurface44P;
  double Asurface48P;
  double Asurface48N;
  double Asurface50P;
  double Asurface52P;
  double Asurface52N;
  double Asurface54P;
  double Asurface54N;
  double Asurface58P;
  double Asurface58N;
  double Asurface60P;
  double Asurface60N;
  double Asurface62P;
  double Asurface62N;
  double Asurface64P;
  double Asurface64N;
  double Asurface88P;
  double Asurface90P;
  double Asurface92P;
  double Asurface92N;
  double Asurface112P;
  double Asurface114P;
  double Asurface116P;
  double Asurface116N;
  double Asurface118P;
  double Asurface118N;
  double Asurface120P;
  double Asurface120N;
  double Asurface122P;
  double Asurface124P;
  double Asurface124N;
  double Asurface206P;
  double Asurface208P;
  double Asurface208N;

  double rsurface0Above;
  double rsurface0Below;
  double rsurface1;
  double asurfaceAbove;
  double asurfaceBelow;
  double EpvolumeAbove;
  double EpvolumeBelow;
  double Avolume0_A;
  double Avolume0_B;
  double Avolume1;
  double Bvolume0Above;
  double Bvolume0Below;
  double Bvolume1;
  double rvolume0Above;
  double rvolume0Below;
  double rvolume1;
  double deltaRvolume;
  double expRvolume;
  double avolume0Above;
  double avolume0Below;
  double avolume1;
  double Ea_above;
  double Ea_below;
  double alphaOverA;
  double Vso;
  double VsoNZ;
  double rso0;
  double rso1;
  double aso0;
  double aso1;
  double AWso;
  double BWso;
  double V_wine;
  double R_wine;
  double rho_wine;
};

struct Parameters {

    Parameters() {}

    double Rc;
    double VHFvol;
    double VHFsur;
    double beta_nl_R0;
    double beta_nl_R1;
    double beta_nl_I0;
    double beta_nl_I1;
    double beta_nl_I0_sur;
    double beta_nl_I1_sur;
    double AHF;
    double RHF;
    double RHFs;
    double aHF;
    double aHFs;
    double fGap_A;
    double fGap_B;
    double BsurfaceA;
    double CsurfaceA;
    double DsurfaceA;
    double Bsurface;
    double Csurface;
    double Dsurface;
    double AsurfaceAbove;
    double AsurfaceBelow;
    double RsurfaceAbove;
    double RsurfaceBelow;
    double EpvolumeAbove;
    double EpvolumeBelow;
    double AvolumeAbove;
    double AvolumeBelow;
    double BvolumeAbove;
    double BvolumeBelow;
    double RvolumeAbove;
    double RvolumeBelow;
    double deltaRvolume;
    double expRvolume;
    double avolumeAbove;
    double avolumeBelow;
    double asurfaceAbove;
    double asurfaceBelow;
    double EaVolume_a;
    double EaVolume_b;
    double alphaVolume;
    double Vso;
    double Rso;
    double aso;
    double AWso;
    double BWso;
    double V_wine;
    double R_wine;
    double rho_wine;
};

Raw_Parameters 
read_par_from_file( const std::string &filename );

Parameters 
get_parameters( const std::string &filename, double A, double Z, double Zp );

NuclearParameters
read_nucleus_parameters( const std::string &filename );

#endif // _READ_PARAMETERS_H_
