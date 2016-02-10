#include "minimizeND.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "reaction.h"



//***************************************************************************
  /**
   *\brief main fitting class
   *
   * performs regular optical model fits to single-energy elastic scattering,
   *
   */
class fitOM: public minimizeND 
{
 public:
  fitOM (string*,string,int,int,bool,bool);
  ~fitOM();
  string title;
  void decodePara();
  void PrintFit();
  void WriteFit(double);
  void SavePara(double *);
  double functND(double*); 
  void grid(int N);
  double gridFit();
  bool outOfBounds();
  static int const TotPara=12; //!< number of  parameters
  int Ndim; //!<Number of actual fit parameters
  int varied[TotPara]; //!<which parameters are fit
  int map1[TotPara];  //!<map from fitted parameters to all parameters
  int map2[TotPara];  //!<map from all parametrs to fitted parameters
  double scaling[TotPara]; //!<scaling of fit parameters to avoid large jumps
  double allPara[TotPara]; //!<values of all the parameters
  bool squared[TotPara]; //!<keep this parameter positive 


  double paraMin[TotPara];
  double paraMax[TotPara];
  double paraBest[TotPara];
  double chiMin;
  int mm;

  // reactions to fit
  int Nreact; //!< number of reactions to fit
  reaction * React[20]; //!< reaction class
  string Rtitle[20]; //!< titles of the reactions

};
