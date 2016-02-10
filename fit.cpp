#include "fit.h"

string const fit::label[TotPara] = 
    { "rC0", "rC1", "VHFvol", "VHFsur", "beta_nl_R0", "beta_nl_R1",
      "beta_nl_I0", "beta_nl_I1", "beta_nl_I0_sur", "beta_nl_I1_sur", "AHF", "rHF0","rHF0s","rHF1","aHF0","aHF0s","aHF1",
      "fGap_A","fGap_B","BsurfaceA","CsurfaceA", "DsurfaceA", "Bsurface","Csurface", "Dsurface", "Asurface9P","Asurface9N",
      "Asurface40P_A","Asurface40N_A", "Asurface40P_B","Asurface40N_B",
      "Asurface42P", "Asurface44P", "Asurface48P", "Asurface48N",
      "Asurface50P", "Asurface52P","Asurface52N", "Asurface54P",
      "Asurface54N", "Asurface58P", "Asurface58N", "Asurface60P",
      "Asurface60N", "Asurface62P", "Asurface62N", "Asurface64P",
      "Asurface64N", "Asurface88P", "Asurface90P", "Asurface92P",
      "Asurface92N", "Asurface112P", "Asurface114P", "Asurface116P",
      "Asurface116N", "Asurface118P", "Asurface118N", "Asurface120P",
      "Asurface120N", "Asurface122P", "Asurface124P", "Asurface124N", 
      "Asurface206P", "Asurface208P","Asurface208N",
      "rsurface0Above","rsurface0Below","rsurface1","asurfaceAbove","asurfaceBelow", "Epvolume_A","Epvolume_B", "Avolume0_A",
      "Avolume0_B", "Avolume_1", "Bvolume0_A","Bvolume0_B", "Bvolume_1", "rvolume0_A","rvolume0_B",
      "rvolume1", "deltaR", "expR", "avolume0_A","avolume0_B", "avolume1", "Ea_Above", "Ea_Below", 
      "alphaOverA", "Vso", "VsoNZ", "rso0","rso1","aso0","aso1","AWso","BWso","V_wine","R_wine","rho_wine"};

//******************************************************************
  /**
   * calculates the total of the chi-squared which is minimized
   \param para is the array of parameter values
   */
double fit::functND(double*para)
{
  
  for (int i=0;i<ND;i++) cout << para[i] << " ";
  cout << endl;
  SavePara(para);

  bool ok;
  for (int i=0;i<Nreact;i++) 
    {

     try
       {
        ok = React[i]->prepare();
       }
     catch(localityException & locExcept)
       {
	 ok = 0;
       } 
     if (ok == 0) return 1e30; //if we cannot adjust fermi level, 
                                //then return large chi
    }
  double out =  0.;
  double chi = 0.;
  for (int i=0;i<Nreact;i++)
    {

     try
       {
       chi=React[i]->ChiSquared();
      }
 
     catch(localityException & locExcept)
       {
         chi = 1000.;
       }

     if (isnan(chi)) chi =1000.;
     out += chi;
    }
/*
 
  //chisq from ratios of total cross sections
  for (int i=0;i<Nratio;i++) 
    {
     double chi = Ratio[i]->ChiSquared();
     if (isnan(chi)) chi=1000.;
     out += chi;
    }
  
*/

  // increase chisq if parameters are outside of reasonable bounds
  out += chiPara();

  cout << "chisq= " << out << endl;

  return out;

}

/**
 *Constructor - reads in experimental data and parameter values
 \param title0 gives name of file with initial parameters
 \param n is number of parameters fit
 */
fit::fit(string *title0,int n) : minimizeND(n)
{
  title = *title0;
  string filename(title + ".inp");
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << endl;
      abort();
    }
  int mvolume, typeN;
  file  >> mvolume;
  file >> typeN;


  string variable;
  for (int i=0;i<TotPara;i++)
    {
      file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;

      if (squared[i]) allPara[i] = sqrt(allPara[i]); 
    }


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





  //read in all the reactions which we will fit
  file.open("reactions.inp");
  file >> Nreact;
  bool yes = 1;
  for (int i=0;i<Nreact;i++)
    {
      file >> Rtitle[i];
      React[i] = new reaction(&(Rtitle[i]),yes,typeN);
      React[i]->mvolume = mvolume;
    }
  file.close();
  file.clear();


 /* 
  //read in all the ratios of total cross section we will fit
  file.open("ratios.inp");
  file >> Nratio;
  for (int i=0;i<Nratio;i++)
    {
      string reaction1;
      string reaction2;
      file >> RatioTitle[i] >> reaction1 >> reaction2;
      int iReact1 = -1;
      int iReact2 = -1;
      for (int j=0;j<Nreact;j++)
	{
          if (reaction1 == Rtitle[j]) iReact1 = j;
          if (reaction2 == Rtitle[j]) iReact2 = j;
	}
      if (iReact1 == -1) cout << "could not find reaction 1" << endl;
      if (iReact2 == -1) cout << "could not find reaction 2" << endl;
      if (iReact1 == -1 || iReact2 == -1) abort();


      Ratio[i] = new ratio(RatioTitle[i],React[iReact1],React[iReact2]);
    }
  file.close();
  file.clear();
*/  

  decodePara();
 
}
//******************************************************************
  /**
   * Destructor
   */
fit::~fit()
{
  for (int i=0;i<Nreact;i++) delete React[i];
}

void fit::decodePara()
{
  /**
   * decodes the parameter file input and send this information to the 
   * classes associated with each potential
   */
  int index = 0;
  rc0 = allPara[index++];
  rc1 = allPara[index++];
  VHFvol = allPara[index++];
  VHFsur = allPara[index++];
  beta_nl_R0 = allPara[index++];
  beta_nl_R1 = allPara[index++];
  beta_nl_I0 = allPara[index++];
  beta_nl_I1 = allPara[index++];
  beta_nl_I0_sur = allPara[index++];
  beta_nl_I1_sur = allPara[index++];
  AHF = allPara[index++];
  rHF0 = allPara[index++];
  rHF0s = allPara[index++];
  rHF1 = allPara[index++];
  aHF0 = allPara[index++];
  aHF0s = allPara[index++];
  aHF1 = allPara[index++];

  fGap_A = allPara[index++];  
  fGap_B = allPara[index++];  
  BsurfaceA = allPara[index++];  
  CsurfaceA = allPara[index++];
  DsurfaceA = allPara[index++];
  Bsurface = allPara[index++];  
  Csurface = allPara[index++];
  Dsurface = allPara[index++];
  Asurface9P = allPara[index++];
  Asurface9N = allPara[index++];
  Asurface40P_A = allPara[index++];
  Asurface40N_A = allPara[index++];

  Asurface40N_A = Asurface40P_A;
  allPara[index-1] = Asurface40N_A;

  Asurface40P_B = allPara[index++];
  Asurface40N_B = allPara[index++];

  Asurface40N_B = Asurface40P_B;
  allPara[index-1] = Asurface40N_B;

  Asurface42P = allPara[index++];
  Asurface44P = allPara[index++];
  Asurface48P = allPara[index++];
  Asurface48N = allPara[index++];
  Asurface50P = allPara[index++];
  Asurface52P = allPara[index++];
  Asurface52N = allPara[index++];
  Asurface54P = allPara[index++];
  Asurface54N = allPara[index++];
  Asurface58P = allPara[index++];
  Asurface58N = allPara[index++];
  Asurface60P = allPara[index++];
  Asurface60N = allPara[index++];
  Asurface62P = allPara[index++];
  Asurface62N = allPara[index++];
  Asurface64P = allPara[index++];
  Asurface64N = allPara[index++];
  Asurface88P = allPara[index++];
  Asurface90P = allPara[index++];
  Asurface92P = allPara[index++];
  Asurface92N = allPara[index++];
  Asurface112P = allPara[index++];
  Asurface114P = allPara[index++];
  Asurface116P = allPara[index++];
  Asurface116N = allPara[index++];
  Asurface118P = allPara[index++];
  Asurface118N = allPara[index++];
  Asurface120P = allPara[index++];
  Asurface120N = allPara[index++];
  Asurface122P = allPara[index++];
  Asurface124P = allPara[index++];
  Asurface124N = allPara[index++];
  Asurface206P = allPara[index++];
  Asurface208P = allPara[index++];
  Asurface208N = allPara[index++];
  rsurface0Above = allPara[index++];
  rsurface0Below = allPara[index++];
  rsurface1 = allPara[index++];
  asurfaceAbove = allPara[index++];
  asurfaceBelow = allPara[index++];
  EpvolumeAbove = pow(allPara[index++],2);
  EpvolumeBelow = pow(allPara[index++],2);
  Avolume0_A = allPara[index++];
  Avolume0_B = allPara[index++];
  Avolume1 = allPara[index++];
  Bvolume0Above = allPara[index++];
  Bvolume0Below = allPara[index++];
  Bvolume1 = allPara[index++];
  rvolume0Above = allPara[index++];
  rvolume0Below = allPara[index++];
  rvolume1 = allPara[index++];
  deltaRvolume = allPara[index++];
  expRvolume = allPara[index++];
  avolume0Above = allPara[index++];
  avolume0Below = allPara[index++];
  avolume1 = allPara[index++];
  Ea_above = allPara[index++];
  Ea_below = allPara[index++];
  alphaOverA = allPara[index++];
  Vso = allPara[index++];
  VsoNZ = allPara[index++];
  rso0 = allPara[index++];
  rso1 = allPara[index++];
  aso0 = allPara[index++];
  aso1 = allPara[index++];
  AWso = allPara[index++];
  BWso = allPara[index++];
  V_wine = allPara[index++];
  R_wine = allPara[index++];
  rho_wine = allPara[index++]; 

  
  for (int i=0;i<Nreact;i++)
    {
     React[i]->Rc = pow(React[i]->A,1./3.)*rc0 + rc1;
     React[i]->VHFvol = VHFvol;
     React[i]->VHFsur = VHFsur;
     //sign = 1 for protons and -1 for neutrons
     //asymmetry = (N-Z)/A

     React[i]->AHF = AHF;
     React[i]->beta_nl_R0 = beta_nl_R0;
     React[i]->beta_nl_R1 = beta_nl_R1;
     React[i]->beta_nl_I0 = beta_nl_I0;
     React[i]->beta_nl_I1 = beta_nl_I1;
     React[i]->beta_nl_I0_sur = beta_nl_I0_sur;
     React[i]->beta_nl_I1_sur = beta_nl_I1_sur;
     React[i]->RHF = pow(React[i]->A,1./3.)*rHF0 + rHF1;
     React[i]->RHFs = pow(React[i]->A,1./3.)*rHF0s+ rHF1;
     React[i]->aHF = aHF0 +aHF1/pow(React[i]->A,1./3.);
     React[i]->aHFs = aHF0s +aHF1/pow(React[i]->A,1./3.);


     React[i]->fGap_A = fGap_A;
     React[i]->fGap_B = fGap_B;
     React[i]->BsurfaceA = BsurfaceA;
     React[i]->CsurfaceA = CsurfaceA;
     React[i]->DsurfaceA = DsurfaceA;
     React[i]->Bsurface = Bsurface;
     React[i]->Csurface = Csurface;
     React[i]->Dsurface = Dsurface;
     if (React[i]->Zp == 1 && React[i]->A == 206)
       {
         React[i]->Asurface_A = Asurface206P;
         React[i]->Asurface_B = Asurface206P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 208)
       {
         React[i]->Asurface_A = Asurface208P;
         React[i]->Asurface_B = Asurface208P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 208)
       {
         React[i]->Asurface_A = Asurface208N;
         React[i]->Asurface_B = Asurface208N;
       }
     else if (React[i]->Zp == 1  && React[i]->A == 40)
       {
         React[i]->Asurface_A = Asurface40P_A;
         React[i]->Asurface_B = Asurface40P_B;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 40)
       {
         React[i]->Asurface_A = Asurface40N_A;
         React[i]->Asurface_B = Asurface40N_B;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 48)
       {
         React[i]->Asurface_A = Asurface48P;
         React[i]->Asurface_B = Asurface48P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 48)
       {
         React[i]->Asurface_A = Asurface48N;
         React[i]->Asurface_B = Asurface48N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 42)
       {
         React[i]->Asurface_A = Asurface42P;
         React[i]->Asurface_B = Asurface42P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 44)
       {
         React[i]->Asurface_A = Asurface44P;
         React[i]->Asurface_B = Asurface44P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 50)
       {
         React[i]->Asurface_A = Asurface50P;
         React[i]->Asurface_B = Asurface50P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 52)
       {
         React[i]->Asurface_A = Asurface52P;
         React[i]->Asurface_B = Asurface52P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 52)
       {
         React[i]->Asurface_A = Asurface52N;
         React[i]->Asurface_B = Asurface52N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 58)
       {
         React[i]->Asurface_A = Asurface58P;
         React[i]->Asurface_B = Asurface58P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 60)
       {
         React[i]->Asurface_A = Asurface60P;
         React[i]->Asurface_B = Asurface60P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 62)
       {
         React[i]->Asurface_A = Asurface62P;
         React[i]->Asurface_B = Asurface62P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 64)
       {
         React[i]->Asurface_A = Asurface64P;
         React[i]->Asurface_B = Asurface64P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 58)
       {
         React[i]->Asurface_A = Asurface58N;
         React[i]->Asurface_B = Asurface58N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 60)
       {
         React[i]->Asurface_A = Asurface60N;
         React[i]->Asurface_B = Asurface60N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 62)
       {
         React[i]->Asurface_A = Asurface62N;
         React[i]->Asurface_B = Asurface62N;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 64)
       {
         React[i]->Asurface_A = Asurface64N;
         React[i]->Asurface_B = Asurface64N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 54)
       {
         React[i]->Asurface_A = Asurface54P;
         React[i]->Asurface_B = Asurface54P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 54)
       {
         React[i]->Asurface_A = Asurface54N;
         React[i]->Asurface_B = Asurface54N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 88)
       {
         React[i]->Asurface_A = Asurface88P;
         React[i]->Asurface_B = Asurface88P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 90)
       {
         React[i]->Asurface_A = Asurface90P;
         React[i]->Asurface_B = Asurface90P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 92)
       {
         React[i]->Asurface_A = Asurface92P;
         React[i]->Asurface_B = Asurface92P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 92)
       {
         React[i]->Asurface_A = Asurface92N;
         React[i]->Asurface_B = Asurface92N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 112)
       {
         React[i]->Asurface_A = Asurface112P;
         React[i]->Asurface_B = Asurface112P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 114)
       {
         React[i]->Asurface_A = Asurface114P;
         React[i]->Asurface_B = Asurface114P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 116)
       {
         React[i]->Asurface_A = Asurface116P;
         React[i]->Asurface_B = Asurface116P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 116)
       {
         React[i]->Asurface_A = Asurface116N;
         React[i]->Asurface_B = Asurface116N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 118)
       {
         React[i]->Asurface_A = Asurface118P;
         React[i]->Asurface_B = Asurface118P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 118)
       {
         React[i]->Asurface_A = Asurface118N;
         React[i]->Asurface_B = Asurface118N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 120)
       {
         React[i]->Asurface_A = Asurface120P;
         React[i]->Asurface_B = Asurface120P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 120)
       {
         React[i]->Asurface_A = Asurface120N;
         React[i]->Asurface_B = Asurface120N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 122)
       {
         React[i]->Asurface_A = Asurface122P;
         React[i]->Asurface_B = Asurface122P;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 124)
       {
         React[i]->Asurface_A = Asurface124P;
         React[i]->Asurface_B = Asurface124P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 124)
       {
         React[i]->Asurface_A = Asurface124N;
         React[i]->Asurface_B = Asurface124N;
       }
     else if (React[i]->Zp == 1 && React[i]->A == 9)
       {
         React[i]->Asurface_A = Asurface9P;
         React[i]->Asurface_B = Asurface9P;
       }
     else if (React[i]->Zp == 0 && React[i]->A == 9)
       {
         React[i]->Asurface_A = Asurface9N;
         React[i]->Asurface_B = Asurface9N;
       }
     else 
       {
         cout << "reaction not implemented " << endl;
         cout << React[i]->Zp << " " << React[i]->A << endl;
	 abort();
       }

     React[i]->RsurfaceAbove = pow(React[i]->A,1./3.)*rsurface0Above + rsurface1;
     React[i]->RsurfaceBelow = pow(React[i]->A,1./3.)*rsurface0Below + rsurface1;


     React[i]->EpvolumeAbove = EpvolumeAbove;
     React[i]->EpvolumeBelow = EpvolumeBelow;

     React[i]->Avolume_A = Avolume0_A + 
       React[i]->asymmetry*React[i]->sign*Avolume1; 

     React[i]->Avolume_B = Avolume0_B + 
       React[i]->asymmetry*React[i]->sign*Avolume1; 

     React[i]->BvolumeAbove = Bvolume0Above + 
       React[i]->asymmetry*React[i]->sign*Bvolume1;
     React[i]->BvolumeBelow = Bvolume0Below + 
       React[i]->asymmetry*React[i]->sign*Bvolume1;



     React[i]->RvolumeAbove = pow(React[i]->A,1./3.)*rvolume0Above + 
           + rvolume1;
     React[i]->RvolumeBelow = pow(React[i]->A,1./3.)*rvolume0Below + 
           + rvolume1;

     React[i]->deltaRvolume = deltaRvolume;

     React[i]->expRvolume = expRvolume;

     React[i]->avolumeAbove = avolume0Above + avolume1/pow(React[i]->A,1./3.);
     React[i]->avolumeBelow = avolume0Below + avolume1/pow(React[i]->A,1./3.);

     React[i]->asurfaceAbove = asurfaceAbove;
     React[i]->asurfaceBelow = asurfaceBelow;
     React[i]->EaVolume_a = Ea_above;
     React[i]->EaVolume_b = Ea_below;
     //alphaOverA = alpha/Avolume
     React[i]->alphaVolume = alphaOverA*React[i]->Avolume_A;

     React[i]->Vso = Vso + React[i]->asymmetry*React[i]->sign*VsoNZ;
     React[i]->Rso = pow(React[i]->A,1./3.)*rso0 + rso1;
     React[i]->aso = aso0 + aso1/pow(React[i]->A,1./3.);
     React[i]->AWso = AWso;
     React[i]->BWso = BWso;
     React[i]->V_wine = V_wine;
     React[i]->R_wine = R_wine;
     React[i]->rho_wine = rho_wine; 

    }
}

//*********************************************************************
  /**
   * saves the fitted and nonfitted parameters in the same form as
   * the input file
   */
void   fit:: SavePara(double * para)
    {
    
  int j_first = map2[0];
  double beta_min = .8;
  double beta_max = 2.0;
  //para[0] = Map_back( beta_min * scaling[j_first], beta_max * scaling[j_first], para[0]);
     for (int i=0;i<Ndim;i++)
       {
         int j = map2[i];
         allPara[j] =para[i]/scaling[j];
        }
  
     decodePara();

    // para[0] = Map_to_inf( beta_min * scaling[j_first], beta_max * scaling[j_first], para[0]);
    
   }
//******************************************************************
  /**
   * Prints to screen the fit parameters
   *
   */
void fit::PrintFit()
{
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

void fit::WriteFit(double chiMin)
{

  //make a linear fit to the AHF values obtained from the main fit
  //fitAHF();

  string filename(title + ".out");
  cout << filename << endl;

  ofstream file(filename.c_str(),ios::out);
  file.setf(ios::fixed); //used fixed precision
  file.precision(11); // set significant digits
  file << setw(10) << React[0]->mvolume << endl;
  file << setw(10) << React[0]->typeN << endl;

  for (int i=0;i<TotPara;i++)
    {
      double value = allPara[i];
      if (squared[i]) value = pow(value,2);  
    file << setw(17) << value << " " << 
            setw(3) << varied[i] << " " << 
            setw(3) << squared[i] << " " <<
            setw(17) << scaling[i] << " " << 
            setw(10) << label[i] << endl;
    }

  file << " " << endl;
  file << " minimun chi squared = " << chiMin << endl;

  file.close();

}
//*************************************************
double fit::chiPara()
{
  double out = 0.;
//  if (asurface < 0.0 || asurface > 8) out += 100.;
//  if (avolume0 < 0.0 || avolume0 > 8) out += 100.;
  if (beta_nl_R0 < 0.0 || beta_nl_R0 > 3) out += 10000.;
  if (beta_nl_R1 < 0.0 || beta_nl_R1 > 3) out += 10000.;
  if (beta_nl_I0 < 0.0 || beta_nl_I0 > 3) out += 10000.;
  if (beta_nl_I1 < 0.0 || beta_nl_I1 > 3) out += 10000.;
  if (asurfaceAbove < 0.0 || asurfaceAbove > 2) out += 10000.;
  if (asurfaceBelow < 0.0 || asurfaceBelow > 2) out += 10000.;
  if (rsurface0Below < 0.0 || rsurface0Below > 3) out += 10000.;
  if (rsurface0Below < 0.0 || rsurface0Below > 3) out += 10000.;
  if (avolume0Above < 0.0 || avolume0Above > 2) out += 10000.;
  if (avolume0Below < 0.0 || avolume0Below > 2) out += 10000.;
  if (rvolume0Above < 0.0 || rvolume0Above > 3) out += 10000.;
  if (rvolume0Below < 0.0 || rvolume0Below > 3) out += 10000.;
  if (aso0 < 0.0 || aso0 > 2) out += 10000.;
  if (rso0 < 0.0 || rso0 > 3) out += 10000.;
  if (Avolume1 < 0.) out += 10000.;

  return out; 
   
}


// Map the parameter p_external in [a, b] to a new parameters
// in [ -Inf, Inf ]
double fit::Map_to_inf( double a, double b, double p_external ) {

    // a is the lower bound
    // b is the upper bound
    return std::asin( 2 * ( p_external - a ) / ( b - a ) - 1 );

}

// Map the parameter p_internal in [-Inf, Inf] to the original parameter 
double fit::Map_back( double a, double b, double p_internal ) {

    // a is the lower bound
    // b is the upper bound
    return a + ( ( b - a ) / 2 ) * ( std::sin( p_internal ) + 1 );
}
