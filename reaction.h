#ifndef _reaction
#define _reaction


#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include "legendre.h"
#include "meshes.h"
#include "numerical.h"
#include "io.h"
#include "Gauss_q.h"


#include "pot.h"
#include "scatterRspace.h"
#include "boundRspace.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>


#include <fstream>
#define root
#ifdef root
  #include "TFile.h"
  #include "TCanvas.h"
  #include "TGraphErrors.h"
  #include "TH2S.h"
  #include "TPad.h"
  #include "TLine.h"
  #include "TLatex.h"
#endif

struct xdata
{
  double energyLab;
  double energyCM;
  double xsec;
  double sigma;
};


struct datas
{
  double energyLab;
  double energyCM;
  string name;
  int nX;
  int nA;
  int fit;
  double* Xtheta;
  double* xsec;
  double* Xsigma;
  double* Atheta;
  double* anal;
  double* Asigma;
};

struct lev
{
  int N;
  double j;
  int l;
  double energy;
};

struct levels
{
  double Energy; // level energy
  double SigmaEnergy;  // experimetal uncertainty
  int N;
  double j;
  int l;
  int color;
  int Efit; // fit level energy
  int Rfit; // fit level rms radius
  int Dfit; // fit level width
  int Sfit; // fit level spectroscopic factor
  double Rrms;     //rms radius
  double SigmaRrms;  //exp uncertainty
  double Delta; // level width (MeV)
  double SigmaDelta; // experimental uncertainty
  double SpectFactor; //spectroscopic Factor
  double SigmaSpect; // error in Spect Factor
};

//calculated level properties
/*
 *\brief stores information on calculated levels
 */
struct levelTh
{
  double N; //!< number of nodes in the radial wavefunction
  double j; //!< total angular momentum of level
  int l; //!< orbital angular momentum of level
  double energy; //!< level energy in MeV
  double Rrms; //!< rms radius of wavefunction
  double SpectFactor; //!< Spectroscopic factor
  double Occupation;  //!<occupation probability of level
  double Width;  //!< width of level in MeV
  int color;
  string name; //!< spectscopic name of level in latex
  double *SpectFunction;
  double SpectE;
  double ANC; //!< asymptotic normalization coefficient
};
int const NdataMax = 40;
int const NlevelMax = 30;


//integrated moment
struct moment
{
  double EnergyLab;
  double EnergyCM;
  double Jreal;
  double Jimag;
  double Jso;
  double RMSreal;
  double RMSimag;
  double RMSso;
};


struct chd_data{
     double r_0;
     double ch_den_0_exp;
     double err_0_exp;
     double r_middle_exp;
     double ch_den_middle_exp;
     double err_middle_exp;
     double R_rms_exp;
     double err_R_rms_exp;
};

  /**
   *\brief deal with a single reaction, ie, p+40Ca
   *
   * performs dispersive optical model calculations of elastic scattering, 
   * reaction and toal sections, bound state energies, RMS radii, widths, and 
   * spectroscopic factors for a single reaction. 
   */


class reaction
{
 public:
  bool DOM; //!< indicates that a DOM fit is taking place
  static double const pi;
  reaction () {};
  reaction (string *title0,bool flag0, int typeN0);
  reaction (string *title0,int jdata,bool btxsec, bool banal);
  ~reaction();
  int Ndata; //!< number of data sets
  double xsecMax;
  double xsecMin;
  int DegreesOfFreedom;
  int LevelDegreesOfFreedom;
  int XsecDegreesOfFreedom;
  bool flag;

  int Nfermi; //!< number of levels defining the Fermi energy
  lev ValenceHole;
  lev ValenceParticle;

  //elastic scattering data
  datas data[NdataMax];
  string title;
  string directory;

  //reaction xsection data
  xdata *Xdata;
  int NXdata;

  //total xsection data
  xdata *TotXdata;
  int NTotXdata;

  bool prepare();
  double ChiSquared(); // returns the chi squared.
  double BoundChiSquared();
  void PlotFit(); // plots data and fitted curve
  void PlotPotentialEcm(); //plots potential as function of Ecm energy 
  void OpenRootFile();
  void CloseRootFile();
  double FindFermi();
  bool AdjustFermi();
  double energyLab2Cm(double);
  double energyCm2Lab(double);
  void Print_DiffXsection( double Elab ); 
  void load();
  void loadOM();

  pot *Pot;
  scatterRspace *ScatterRspace;
  boundRspace * BoundRspace;


  double Zp; //charge of projectile 
  double Z;
  double A;
  double asymmetry; // (N-Z)/A
  double sign;
  double Rc;
  double Efermi;
  double gap; //!< gap of particle sepcies under consideration
  double gapOther; //!< gap for other particle species (proton or neutron)
  double gapMin; //!< minimum of proton and neutron gaps
  double Wstart; //!< energy above and below Efermi where Imaginary pot starts
  double EpvolumeAbove;
  double EpvolumeBelow;
  int readCoulomb; //!< specify whether or not to read in exp. charge density

  int mvolume;

  int Ndim;

  double VHFvol;
  double VHFsur;
  double AHF;
  double beta_nl_R0; //!< nonlocality distance for real potential
  double beta_nl_R1; //!< nonlocality distance for real potential
  double beta_nl_I0; //!< nonlocality distance for imag potential
  double beta_nl_I1; //!< nonlocality distance for imag potential
  double beta_nl_I0_sur; //!< nonlocality distance for imag potential
  double beta_nl_I1_sur; //!< nonlocality distance for imag potential
  double RHF;
  double aHF;
  double RHFs;
  double aHFs;


  double fGap_A;
  double fGap_B;
  double Asurface_A; // Surface strength above Fermi energy
  double Asurface_B; // Surface strength below Fermi energy
  double BsurfaceA;
  double CsurfaceA;
  double DsurfaceA;
  double Bsurface;
  double Csurface;
  double Dsurface;

  double RsurfaceAbove;
  double RsurfaceBelow;
  double asurfaceAbove;
  double asurfaceBelow;

 
  double Avolume_A; // Volume strength above Fermi energy
  double Avolume_B; // Volume sterngth below Fermi energy
  double BvolumeAbove;
  double BvolumeBelow;
  double deltaRvolume;
  double expRvolume;
  double RvolumeAbove;
  double RvolumeBelow;
  double evolume;
  double avolumeAbove;
  double avolumeBelow;
  int AsyVolume;
  double alphaVolume;
  double EaVolume_a;
  double EaVolume_b;


  double Vso;
  double Rso;
  double aso;
  double AWso;
  double BWso;
  int typeN;
  double V_wine; 
  double R_wine;
  double rho_wine;
 
  levels Level[NlevelMax];
  int Nlevel;
  int NFitLevels; // number of levels in fit

  levelTh LevelTh[NlevelMax];
  levelTh LevelThSort[NlevelMax];
  int NlevelTh;
  double Ef; // calculated Fermi energy
  moment Moment[100];
  int Nmoment;

  
  chd_data Chi_Squared_data;
  // spectrospopic function
  int Nsf; // number of energy in array
  double Elow; // for electon scattering lower limits for integrated strength
  double Ehigh; // upper limit
  double *Esf; //array of energies

  double chidif;
  ofstream Fout;
#ifdef root
  TFile *f;
#endif


 private:
  bool bprint;
 
};
#endif
