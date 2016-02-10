#include "minimizeND.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "reaction.h"
//#include "ratio.h"
#include "localityException.h"

//***************************************************************************
  /**
   *\brief main fitting class
   *
   * performs dispersive optical model fits to elastic scattering, reaction 
   * and toal sections, bound state energies, RMS radii, widths, and 
   * spectroscopic factors
   */
class fit: public minimizeND 
{
 public:
  fit (string*,int);
  ~fit();
  string title;
  void decodePara();
  void PrintFit();
  void WriteFit(double);
  void SavePara(double *);
  double functND(double*); 


  static int const TotPara=101; //!< number of  parameters
  int Ndim; //!<Number of actual fit parameters
  int varied[TotPara]; //!<which parameters are fit
  int map1[TotPara];  //!<map from fitted parameters to all parameters
  int map2[TotPara];  //!<map from all parametrs to fitted parameters
  double scaling[TotPara]; //!<scaling of fit parameters to avoid large jumps
  double allPara[TotPara]; //!<values of all the parameters
  bool squared[TotPara]; //!<keep this parameter positive 
  void fitAHF(); 
  double chiPara();
  double Map_to_inf( double, double , double);
  double Map_back(double , double , double);
  // reactions to fit
  int Nreact; //!< number of reactions to fit
  reaction * React[30]; //!< reaction class
  string Rtitle[30]; //!< titles of the reactions

  // ratios of total cross section to fit
  int Nratio; //!< number of ratios
//  ratio * Ratio[2]; //!< ratio class
//  string RatioTitle[2]; //!< titles of the ratios

 private:
  static string const label[TotPara]; 



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
  double Asurface206P;
  double Asurface208P;
  double Asurface208N;
  double Asurface40P_A;
  double Asurface40N_A;
  double Asurface40P_B;
  double Asurface40N_B;
  double Asurface48P;
  double Asurface48N;
  double Asurface42P;
  double Asurface44P;
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
  double Asurface9P;
  double Asurface9N;
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
  int typeN;
};
