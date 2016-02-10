#include "fitOM.h"

//******************************************************************
  /**
   * calculates the total of the chi-squared which is minimized
   *\para para is the array of parameter values
   */
double fitOM::functND(double*para)
{
  //for (int i=0;i<ND;i++) cout << para[i] << " " ;
  //cout << endl;
  SavePara(para);
  React[0]->loadOM();
  double out =  0.;
  out = React[0]->ChiSquared();
  //cout << "chisq= " << out << endl;

  //check if a parameter is out of the bounds listed in the *.inp file
  if (outOfBounds()) out += 100.;
  
  if (isnan(out))
    {
     cout << "chisq is nan" << endl;
     out = 1e20;
     //abort();
    }
  return out;
}

/**
 *Constructor - reads in experimental data and parameter values
 \param title0 gives name of file with initial parameters
 \param reactionTitle is name of reaction
 \param we fir the (jdata)th data set in file reactionTitle.data  
 \para n is number of parameters fit
 */
fitOM::fitOM(string *title0, 
	     string reactionTitle, int jdata, int n, bool btxsec, bool banal) 
             : minimizeND(n)
{


  Nreact = 1;
  React[0] = new reaction(&reactionTitle,jdata,btxsec,banal);
  if (React[0]->data[0].nX <= 0) return;

  title = *title0;
  string filename(title + ".inp");
  cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << endl;
      abort();
    }

  string variable;
  mm = 0;
  for (int i=0;i<TotPara;i++)
    {
      file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] 
           >> paraMin[i] >> paraMax[i] >> variable;
      if (i == 4 && React[0]->data[0].energyCM-React[0]->Efermi > 50.)
	{ 
          varied[i] = 0;
          allPara[i] = 0.;
	}
      if (i == 5 && React[0]->data[0].energyCM-React[0]->Efermi < 25.)
	{ 
          varied[i] = 0;
          allPara[i] = 0.;
	}
      if (squared[i]) allPara[i] = sqrt(allPara[i]); 
      if (varied[i] == 1) mm++;
    }
  if (mm != n) return; 

  //create maps of all paremeter to fitted parameters
  Ndim = 0;
  for (int i=0;i<TotPara;i++)
    {
      map1[i] = -1;
      map2[i] = -1;
      if (varied[i]) 
	{
	 map1[i] = Ndim;
         map2[Ndim] = i;
         Ndim++;
	}
    }
  file.close();
  file.clear();





  decodePara();
 
}
//******************************************************************
  /**
   * Destructor
   */
fitOM::~fitOM()
{
  delete React[0];
}

void fitOM::decodePara()
{
  /**
   * decodes the parameter file input and send this information to the 
   * classes associated with each potential
   */
  int index = 0;
  double rc = allPara[index++];
  double VHF = allPara[index++];
  double rHF = allPara[index++];
  double aHF = allPara[index++];

  double Asurface = allPara[index++];
  double Avolume = allPara[index++];
  double rimag = allPara[index++];
  double aimag = allPara[index++];

  double Vso = allPara[index++];
  double Wso = allPara[index++];
  double rso = allPara[index++];
  double aso = allPara[index++];

  React[0]->Rc = pow(React[0]->A,1./3.)*rc;
  React[0]->VHFvol = VHF;
  React[0]->RHF = pow(React[0]->A,1./3.)*rHF;
  React[0]->aHF = aHF;

  React[0]->Asurface = Asurface;
  React[0]->Rsurface = pow(React[0]->A,1./3.)*rimag;
  React[0]->asurface = aimag;

  React[0]->Avolume = Avolume; 
  React[0]->Rvolume = pow(React[0]->A,1./3.)*rimag;
  React[0]->avolume = aimag;

  React[0]->Vso = Vso;
  React[0]->Rso = pow(React[0]->A,1./3.)*rso;
  React[0]->aso = aso;
  React[0]->AWso = Wso;
}

//*********************************************************************
  /**
   * saves the fitted and nonfitted parameters in the same form as
   * the input file
   */
void   fitOM:: SavePara(double * para)
    {
     for (int i=0;i<Ndim;i++)
       {
         int j = map2[i];
         allPara[j] = para[i]/scaling[j];
        }

     decodePara();
    }
//******************************************************************
  /**
   * Prints to screen the fit parameters
   *
   */
void fitOM::PrintFit()
{
  string label[TotPara] = {"rC","VH","rHF","aHF","Asurface",
                           "Avolume","rimag","aimag",
                           "Vso","Wso",
                           "rso","aso"};
  for (int i=0;i<TotPara;i++)
    {
      double value = allPara[i];
      if (squared[i]) value = pow(value,2);  
      cout << setw(17) << value << " " << 
             setw(3) << varied[i] << " " << 
             setw(3) << squared[i] << " " <<
             setw(17) << scaling[i] << " " << 
             setw(10) << label[i] << endl;
    }
}
//******************************************************************
  /**
   * writes out the fitted parameter in the file title.out
   * This file can later renamed title.inp and used as the 
   * input file for another fit 
   *
   */

void fitOM::WriteFit(double chiMin)
{


  string filename(title + ".out");
  cout << filename << endl;
  ofstream file(filename.c_str(),ios::out);
  file.setf(ios::fixed); //used fixed precision
  file.precision(11); // set significant digits
  //file << setw(10) << React[0]->msurface << " " << 
  //        setw(10) << React[0]->mvolume << endl;

  string label[TotPara] = {"rC","VH","rHF","aHF","Asurface",
                           "Avolume","rimag","aimag",
                           "Vso","Wso",
                           "rso","aso"};


  for (int i=0;i<TotPara;i++)
    {
      double value = allPara[i];
      if (squared[i]) value = pow(value,2);  
    file << setw(17) << value << " " << 
            setw(3) << varied[i] << " " << 
            setw(3) << squared[i] << " " <<
            setw(17) << scaling[i] << " " <<
            setw(17) << paraMin[i] << " " <<
            setw(17) << paraMax[i] << " " <<
            setw(10) << label[i] << endl;
    }

  file << " " << endl;
  file << " minimun chi squared = " << chiMin << endl;

  file.close();

}


//conducts a grid search

void fitOM::grid(int N)
 {
   int const Nsteps =5;
   for (int i=0;i<Nsteps;i++)
     {
       if (varied[N] && N < 3) cout << N << " " << i << endl;

       if (varied[N])
	 allPara[N] = paraMin[N] + (double)i/(double)(Nsteps-1)*(paraMax[N]-paraMin[N]);
      
      if (N+1 == TotPara)
        {
          decodePara();
          React[0]->loadOM();
          double chi =  0.;
          chi = React[0]->ChiSquared();
          //cout << "chisq= " << chi << " " << chiMin << endl;
          if (isnan(chi))
	    {
              //abort();
              chi = 1.e20;
	    }
          if (chi < chiMin)
            {
              chiMin = chi;
              cout << chiMin << endl;
              for (int j=0;j<TotPara;j++) paraBest[j] = allPara[j];
            }
        }
      else grid(N+1);
      if (!varied[N]) break;
     }
 }

//**************************************************************
double fitOM::gridFit()
{
  chiMin = 1000000000000.;
  grid(0);
  for (int i=0;i<TotPara;i++) allPara[i] = paraBest[i];
  return chiMin;
}
//***************************************************************
bool fitOM::outOfBounds()
{
  for (int i=0;i<TotPara;i++)
    {
      if (allPara[i] < paraMin[i]) return 1;
      if (allPara[i] > paraMax[i]) return 1;
    }

  return 0;
}
