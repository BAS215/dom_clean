// to run this program 
// chisq X  to fit data
// chisq a  to not fit

#include "fit.h"


int main(int Narg, char* pszArg[])
{

  if (Narg > 3) return 0;
  string title(pszArg[1]);
  string fileName = title+".inp";

  //read in input file to cound how many fitted parameters
  ifstream file(fileName.c_str());
  if (file.fail())
    {
      cout << "could not open file " << fileName << endl;
      abort();
    }
  char line[80];
  //skip two lines
  file.getline(line,80);
  file.getline(line,80);
  string variable;
  double parameter,varied,squared,scaling;
  int Ndim = 0;
  for (int i=0;i<fit::TotPara;i++)
    {
      file >> parameter >> varied >> squared>>scaling >> variable;
      if (varied) Ndim++;
    }  
  

  file.close();
  file.clear();




  fit Fit(&title,Ndim);
  
  double * para;
  para = new double [Ndim];

  for (int i=0;i<Ndim;i++)
    {
      int j= Fit.map2[i];
      para[i] = Fit.allPara[j]*Fit.scaling[j];
    }
  int j_first=Fit.map2[0];
  double beta_min=.8;
  double beta_max=2.0;
  //para[0]=Fit.Map_to_inf( beta_min * Fit.scaling[j_first], beta_max * Fit.scaling[j_first], para[0]);

  double* xi[Ndim];
  for (int i=0;i<Fit.ND;i++) 
    {
      xi[i] = new double [Ndim];
      for (int j=0;j<Fit.ND;j++)
	{
	  if (i == j) xi[i][j] = 1.;
	  else xi[i][j] = 0.;
	}
    }
  double const ftol = 1.e-2;
  double chiMin=0.;

  string toFit(pszArg[2]);  
  //string toFit("X");
  cout << toFit << endl;
  if (toFit == "X" || toFit ==  "x")chiMin=Fit.powell(para,xi,ftol);

  cout << "minimum chisq= " << chiMin << endl;

  Fit.SavePara(para);

  Fit.WriteFit(chiMin);

  for (int i=0;i<Fit.Nreact;i++)
    {
     Fit.React[i]->OpenRootFile();
     cout << "PlotFit" << endl;
     Fit.React[i]->PlotFit();
     cout << "plotPot" << endl;
     Fit.React[i]->PlotPotentialEcm();

     double Elab = 26;
     Fit.React[i]->Print_DiffXsection( Elab );

      	

     /*
     //write  out transmission coefficients
     fileName = Fit.Rtitle[i]+".CN";
     file.open(fileName.c_str());
     file >> fileName;
     file >> fileName;
     Fit.React[i]->PrintTransCoef(fileName);
     file.close();
     file.clear();
     */

     Fit.React[i]->CloseRootFile();

    }
 

  /*
  for (int i=0;i<Fit.Nratio;i++)
    {
     Fit.Ratio[i]->OpenRootFile();
     cout << "PlotFit" << endl;
     Fit.Ratio[i]->PlotFit();
     Fit.Ratio[i]->CloseRootFile();

    }
    */ 
  delete [] para;
  for (int i=0;i<Ndim;i++)
    {
     delete [] xi[i];
    }
  return 0;
}

