#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "reaction.h"


#define root

#ifdef root
#include "TH2S.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#endif

using namespace std;

/**
 * structure to store experimental ratios of xcestions
 */

struct rdata
{
  double energyLab;  //!< lab energy in MeV
  double energyCM1;  //!< center-of-mass energy of reaction 1 in MeV
  double energyCM2; //!< center-of-mass energy of reaction 2 in MeV
  double x;        //!< experimental ratio
  double sigma;    //!< stand error for experimental data
};


/**
 *\brief deals with fitting ratios of cross sections
 *
 *This class calculates the ratio of total cross sections for two
 *neutron induced reactions.\f$ \frac{\sigma_2-\sigma_1}{\sigma_2+\sigma_1} \f$
 * It also calculates the chi squared relative
 *to experimental data read in.
 */

class ratio
{
 protected:
  string title;  //!< root title of input files

  //ratio data
  rdata *Rdata;  //!< storage of experimental ratios
  int Ndata;     //!< number of experimental data
  reaction * Reaction1;  //!< first reaction for ratio
  reaction * Reaction2;  //!<second reaction for ratio


 public:
  ratio(string title0, reaction *Reaction10, reaction* Reaction20);
  ~ratio();
  double ChiSquared();
  double getRatio(double Elab);
  void PlotFit();
  void OpenRootFile();
  void CloseRootFile();

#ifdef root
  TFile *f;    //!< root file class
#endif


};
