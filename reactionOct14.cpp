#include "reaction.h"

double const reaction::pi = acos(-1.);
bool nonloc2 = true;

//*********************************************************************
  /**
   *constructor reads in data from file
   * Set falg to zero is class used only to gives optical model potential
   \param title0 for input data
   \param flag0 indicates integrating of wavefunctions takes place
   \param typeN0 - type of nonlocality 
  */
reaction::reaction(string *title0,bool flag0, int typeN0)
{
  typeN = typeN0;
  bprint = 1;
  for (int i=0;i<NlevelMax;i++) LevelTh[i].SpectFunction = NULL;
  Esf = NULL;
  NXdata = 0;
  NTotXdata = 0;
  DOM = 1;
  flag = flag0;
  //directory = "/home/bob/DOMA/";
  directory = "Data/";
  string input_dir = "Input/";


  string line;
  title = *title0;
  string filename(input_dir + title + ".inp");
  if (flag) cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << " " << title << endl;
      abort();
    }

  file >> Zp >> Z >> A >> Efermi >> readCoulomb;
  file >> Nfermi;
  if (Nfermi < 1 || Nfermi > 2) cout << "Nfermi not possible" << endl;
  file >> ValenceHole.N >> ValenceHole.j >> ValenceHole.l >> 
          ValenceHole.energy;
  if (Nfermi == 2) 
     file >> ValenceParticle.N >> ValenceParticle.j >> ValenceParticle.l
          >> ValenceParticle.energy;
  file >> gap  >> gapOther;
  gapMin = min(gap,gapOther);
  Wstart = gap/2. + gapMin;

  //Wstart = 0.9*gap;


  asymmetry = (A-2.*Z)/A;
  if (Zp == 1.) sign = 1;
  else if (Zp == 0.) sign = -1.;
  else cout << " unknown Zp " << endl;
  if (file.bad()) cout << "problem with input file" << endl;
  if (file.eof()) cout << "eof with input file" << endl;

  file.close();
  file.clear();


  //asymmetry
  AsyVolume = 1;  // asymmetry for volume
  //alphaVolume = 1.65;
  //EaVolume = 140.;
  //EaVolume = 60.;
  double rmax = 12;
  int ham_pts = 180; // number of points used in construction of hamiltonian
  int Lmax = 5;
  
  //create class for nonlocal potential energy
  //int typeN = 1; // set to U((r1+r2)/2) form
  //int typeN = 0; // set to U*sqrt(f(r1)*f(r2)) form
   int typeN = 1;
  Pot = new pot(typeN);
  Pot->init(Z, Zp, A, readCoulomb ); //initialize
  prepare(); //load potential parameters into Pot

  // constructor for Rspace scttering calculations
  ScatterRspace = new scatterRspace(Z,Zp,A,Pot,title0);
  BoundRspace = new boundRspace( rmax, ham_pts, Efermi, gap, Lmax, Z,Zp,A,Pot);
 
  
  if (flag == 0) return;

  //___________________________________________________________________
  //read elastic scattering data


  int iData = 0;
  Ndata = 0;
  DegreesOfFreedom = 0;
  xsecMax = 0.;
  xsecMin = 1.E32;

  //open data file
  filename = directory + title + ".data";
  cout << filename << endl;
  ifstream fileData (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileData.fail()) 
    {
      cout << "couldn't open data file" << fileData << endl;
    }

  else
    {
      for(;;)
	{
	  int n;
	  string xa;
	  double Elab;
	  int nfit;
	  fileData >> Elab >> nfit;
	  // if nfit ==1 data is included in fit, =0 plotted but not fitted
	  //=-1 data is ignored
          //=2  data is ratio to ruth and in fit
	  if (Elab < 0.) break;


	  getline(fileData,line);
	  getline(fileData,line);
	  cout << line << endl;

	  if (nfit >= 0)
	    {
	      data[iData].energyLab = Elab;
              data[iData].name = line;
	      data[iData].fit = nfit;
	      double Ecm = energyLab2Cm(Elab);
	      data[iData].energyCM = Ecm;
              // if we need to calculate rutherford, load in the Ecm to 
	      ScatterRspace->newEnergy(Ecm,Elab);
	    }

	  if (fileData.eof()) break;
	  if (fileData.bad())
	    {
	      cout << "trouble reading file" << endl;
	      break;
	    }

	  if (nfit >= 0) Ndata++;

	  //make room for cross xsection data
	  fileData >> n >> xa;
          //cout << "n = " << n << " xa = " << xa << endl;
	  if (xa != "X") cout << "X  problem for Elab = " << Elab << endl;

	  if (nfit >= 0) data[iData].nX = n;
	  if (n > 0)
	    {  

	      if (nfit >= 0)
		{
		  data[iData].Xtheta = new double[n];
		  data[iData].xsec = new double[n];
		  data[iData].Xsigma = new double[n];
		}

	      double theta, xsec, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> xsec >> sigma;
     
		  if (nfit >= 0)
		    {

		      data[iData].Xtheta[i] = theta;
		      data[iData].xsec[i] = xsec;
		      data[iData].Xsigma[i] = sigma;
                      


                      // input data is given as ratio to rutherford
                      if( nfit == 2) 
			{
                         double ruth =ScatterRspace->Rutherford(theta*pi/180.); 
                         data[iData].xsec[i]*= ruth;
			
			 if (sigma > 0) data[iData].Xsigma[i] *= ruth;
			}		    


		      if (data[iData].Xsigma[i] < 0.) 
			data[iData].Xsigma[i] *= -data[iData].xsec[i]/100.;

		      if (data[iData].Xsigma[i] < data[iData].xsec[i]*.1)
			data[iData].Xsigma[i] = data[iData].xsec[i]*.1;

                      if ( Zp == 0 && A == 92 && Elab > 10 ) 
                      data[iData].Xsigma[i] /= 5.;
                      if ( Zp == 1 && Elab > 14. && Elab < 100 && 
			   (A==48 || A==44))
                      data[iData].Xsigma[i] /= 5.;

                      if ( Elab > 14. && Elab < 50 && Z == 28)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 1 && Elab > 100. && A==9)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 1 && Elab > 40. && Elab < 100. && A==208)
                      data[iData].Xsigma[i] /= 5.;


                      //if ( Zp == 1 && Elab > 40. && Elab < 100. && A==92)
                      //data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 0 && Elab > 20. && Elab < 40. && A==208)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 0 && Elab > 20. && Elab < 40. && A==9)
                      data[iData].Xsigma[i] /= 5.;

		      if (nfit) DegreesOfFreedom++;

		      xsecMax = max(data[iData].xsec[i],xsecMax);
		      xsecMin = min(data[iData].xsec[i],xsecMin);
		    }
		}
	    }
	  //analysing power data
	  fileData >> n >> xa;

	  if (xa != "A") cout << "A problem for Elab = " << Elab << endl;
	  if (nfit >= 0) data[iData].nA = n;
	  if (n > 0)
	    {
	      if (nfit >= 0)
		{
		  data[iData].Atheta = new double[n];
		  data[iData].anal = new double[n];
		  data[iData].Asigma = new double[n];
		}
	      double theta, pol, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> pol >> sigma;

		  if (nfit >=0)
		    {
		      data[iData].Atheta[i] = theta;
		      data[iData].anal[i] = pol;
		      if (sigma < 0.05) sigma = 0.05; //****** ******
		      data[iData].Asigma[i] = sigma;
                      if (data[iData].energyLab > 50 && data[iData].energyLab
			  < 100) data[iData].Asigma[i] /= 4.;
		      DegreesOfFreedom++;
		    }
		}
	    }

 
	  if (nfit >=0 )
	    {
	      iData++;
	      if (Ndata > NdataMax) 
		{
		  cout << "increase NdataMax" << endl;
                  abort();
		  break;
		}
	    }
	} 

      cout << Ndata << " blocks of data read from file " << filename << endl;
      cout << " degrees of freedom= " << DegreesOfFreedom << endl;
      cout << " Max xsection = " << xsecMax << endl;
      cout << " Min xsection = " << xsecMin << endl;

      fileData.close();
    }
//////////////////////////////////////////////
//////////////////////////////////////////////

  //--------------------- read in levels
  //read in levels
  filename = directory + title + ".lev";
  cout << filename << endl;
  ifstream fileLevel (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileLevel.fail()) 
    {
      cout << "couldn't open data file" << endl;
      NFitLevels = 0;
      Nlevel = 0;
     LevelDegreesOfFreedom = 0;
    }
  else
    {
      getline(fileLevel,line);
     cout << line << endl;
      getline(fileLevel,line);
     //cout << line << endl;
     int Ilevel=0;

     LevelDegreesOfFreedom = 0;

     NFitLevels = 0;
     for (;;)
       {
           fileLevel >> Level[Ilevel].Energy >> Level[Ilevel].SigmaEnergy >>
           Level[Ilevel].N >>
	   Level[Ilevel].j >> Level[Ilevel].l >> Level[Ilevel].color >>
   	   Level[Ilevel].Efit >>
           Level[Ilevel].Rrms >> Level[Ilevel].SigmaRrms >>
	   Level[Ilevel].Rfit >>
           Level[Ilevel].Delta >> Level[Ilevel].SigmaDelta >> 
	   Level[Ilevel].Dfit >>
           Level[Ilevel].SpectFactor >> Level[Ilevel].SigmaSpect >>
           Level[Ilevel].Sfit;

	   //if (A >= 40.) Level[Ilevel].SigmaSpect /= 2.;
        

         if (fileLevel.eof()) break;


         if (Level[Ilevel].Efit) 
	   {
            LevelDegreesOfFreedom++;
	    NFitLevels++;
	   }
         if (Level[Ilevel].Rfit) LevelDegreesOfFreedom++;
         if (Level[Ilevel].Dfit) LevelDegreesOfFreedom++;
         if (Level[Ilevel].Sfit) LevelDegreesOfFreedom++;


         if (Ilevel == NlevelMax-1) 
	   {
	     cout << "increase NlevelMax" << endl;
             Ilevel++;
	     break;
	   } 
         Ilevel++;
       }
     Nlevel = Ilevel;
     cout << Nlevel << " levels read in " << endl;
     fileLevel.close();
    }



  //prepare for spectral functions
  Nsf = 440;
  for (int i=0;i<NlevelMax;i++) LevelTh[i].SpectFunction = new double [Nsf];
  Elow = Efermi - 22.;
  Ehigh = Efermi;
  Esf = new double[Nsf];
  for (int i=0;i<Nsf;i++) Esf[i] = -(double)i/4.+30;



  XsecDegreesOfFreedom = 0;

  //readin reaction xsection data --------
 filename = directory+title + ".xsec";
  cout << filename << endl;
  file.open(filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open reaction *.xsec file" << endl;
      cout << " no reaction xsec in fit" << endl;
      NXdata = 0;
    }
  else 
    {
      getline(file,line);
      getline(file,line);
     file >> NXdata;
     XsecDegreesOfFreedom += NXdata;
     Xdata = new xdata [NXdata];

     for (int i=0;i<NXdata;i++)
       {
         file >> Xdata[i].energyLab >> Xdata[i].xsec >> Xdata[i].sigma;
         double Elab = Xdata[i].energyLab;
         double Ecm = energyLab2Cm(Elab);
         Xdata[i].energyCM = Ecm;
       }
    }
  file.close();
  file.clear();

  cout << "xsec points = " << NXdata << endl;
  if (Zp == 0.) //for neutrons total cross section data
    {
      filename = directory+title + ".txsec";
      cout << filename << endl;
      file.open(filename.c_str(),ios::in);
      // if one cannot open file quit
      if (file.fail()) 
	{
	  cout << "couldn't open reaction *.txsec file" << endl;
	  cout << " no total xsec in fit" << endl;
	  NTotXdata = 0;
	}
      else 
	{
          getline(file,line);
          getline(file,line);

	  file >> NTotXdata;
	  XsecDegreesOfFreedom += NTotXdata;
	  TotXdata = new xdata [NTotXdata];

	  for (int i=0;i<NTotXdata;i++)
	    {
	      file >> TotXdata[i].energyLab >> 
              TotXdata[i].xsec >> TotXdata[i].sigma;

	      //make sure the error bars are not too small
              if (TotXdata[i].sigma < TotXdata[i].xsec*.02)
		TotXdata[i].sigma = 0.02*TotXdata[i].xsec;


             
	      double Elab = TotXdata[i].energyLab;
	      double Ecm = energyLab2Cm(Elab);
              //if (A == 208 && Elab > 25 ) TotXdata[i].sigma /= 5;
	      TotXdata[i].energyCM = Ecm;
              if (TotXdata[i].sigma < 0.) TotXdata[i].sigma *= 
					    -TotXdata[i].xsec/100.;
	    }
	}
    }
  else NTotXdata = 0;
  cout << "tot xsec points" << NTotXdata << endl;

  file.close();
  file.clear();
  //readin integrated moments
 filename = directory+title + ".Vint";
  cout << filename << endl;
  file.open(filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open integtaed moment *.Vint file" << endl;
      cout << " no moments in fit" << endl;
      Nmoment = 0;
    }
  else 
    {
       getline(file,line);
      int flag;
      file >> Nmoment >> flag;
      cout << Nmoment << " " << flag << endl;
      if (flag == 0) Nmoment = 0;

      if (Nmoment > 0)
	{
        string ref;
	for (int i=0;i<Nmoment;i++)
	  {
	    file >> Moment[i].EnergyLab >> Moment[i].Jreal >>
	      Moment[i].Jimag >> Moment[i].RMSreal >> Moment[i].RMSimag >>
	      Moment[i].Jso >> Moment[i].RMSso;
	    getline(file,ref);
            Moment[i].EnergyCM = energyLab2Cm(Moment[i].EnergyLab);

	    cout << Moment[i].EnergyLab << " " << Moment[i].EnergyCM << endl;
	  }
	}
    }
  file.close();
  file.clear();

    // read in charge distribution

//reading boundrspace data to calculate chi squared//
std::ifstream filein("ca40_bound_chd.dat");

filein >> Chi_Squared_data.r_0 >> Chi_Squared_data.ch_den_0_exp 
       >> Chi_Squared_data.err_0_exp >> Chi_Squared_data.r_middle_exp
       >> Chi_Squared_data.ch_den_middle_exp >> Chi_Squared_data.err_middle_exp
       >> Chi_Squared_data.R_rms_exp >> Chi_Squared_data.err_R_rms_exp;

  filein.close();
}

//******************************************************************
  /**
   * Destructor
   */
reaction::~reaction()
{
  delete BoundRspace;
  delete ScatterRspace;
  delete Pot;

  cout << "destroying reaction" << endl;
  if (flag == 0) return;
  for (int i=0;i<Ndata;i++)
    {
      if (data[i].nX > 0)
	{
         delete [] data[i].Xtheta;
         delete [] data[i].xsec;
         delete [] data[i].Xsigma;
	}
      if(data[i].nA > 0)
	{
          delete [] data[i].Atheta;
	  delete [] data[i].anal;
	  delete [] data[i].Asigma;
	}
    }
  for (int i=0;i<NlevelMax;i++) 
    {
      if (LevelTh[i].SpectFunction == NULL) continue;
     delete  [] LevelTh[i].SpectFunction;
    }
  delete [] Esf;
  if (NXdata > 0) delete [] Xdata;
  if (NTotXdata > 0) delete [] TotXdata;
  cout << " reaction distroyed" << endl;
}
//*****************************************************************
  /**
   * loads in in new parameters for the potential
   */
bool reaction::prepare()
{
  Pot->load( Rc, VHFvol, VHFsur, RHF, aHF,  RHFs, aHFs, beta_nl_R0, AHF,
             beta_nl_R1, RsurfaceAbove, RsurfaceBelow, asurfaceAbove, asurfaceBelow, Asurface_A,
             Asurface_B, BsurfaceA, CsurfaceA, DsurfaceA, Bsurface, Csurface, Dsurface,
	     Wstart*fGap_A, Wstart*fGap_B, Efermi, beta_nl_I0, beta_nl_I1, beta_nl_I0_sur, beta_nl_I1_sur,
             RvolumeAbove,RvolumeBelow, deltaRvolume, expRvolume,
             avolumeAbove,avolumeBelow, Avolume_A, Avolume_B, BvolumeAbove, BvolumeBelow ,
             EpvolumeAbove,EpvolumeBelow, mvolume, AsyVolume, alphaVolume, EaVolume_a,
             EaVolume_b, Rso, aso, Vso, AWso, BWso , V_wine , R_wine , rho_wine);

  return true;
  return AdjustFermi();
}

//******************************************************************
  /**
   * returns chi squared per degree of freedom of the fit
   */


double reaction::ChiSquared()
{
  double sum = 0.;
  double sumXsec = 0.;
  //loop over number of energies in fit
  for(int i=0;i<Ndata;i++)
    {
      if (data[i].fit == 0) continue;
      double Ecm =data[i].energyCM;
      double Elab =data[i].energyLab;

      // integrate wavefunction and find phaseshifts
      if (ScatterRspace->getSmatrix(Ecm,Elab)==0) return 1.e6;
       //add prepare compound elastic
       if (Ecm < 15.) ScatterRspace->statistical(ScatterRspace->konst,Ecm);

       // loop over angles
       for (int j=0;j<data[i].nX;j++)
	 {
	   double angle = data[i].Xtheta[j]*pi/180.;
           double xsecTheory = ScatterRspace->DifferentialXsection(angle);
           // add in compound elastic
           if (Ecm < 15.) xsecTheory += ScatterRspace->DifferentialXsectionCE(angle);
           //add to chi squared
           sum += pow(xsecTheory-data[i].xsec[j],2)/
	     pow(data[i].Xsigma[j],2);

	 }
       for (int j=0;j<data[i].nA;j++)
	 {
	   double angle = data[i].Atheta[j]*pi/180.;
           double xsecTheory = ScatterRspace->DifferentialXsection(angle);
           //add to chi squared

           double A = ScatterRspace->AnalyzePower;
	   if (Ecm < 15.) 
	     {
               double shape = ScatterRspace->DifferentialXsectionCE(angle);
               A *= xsecTheory/(xsecTheory+shape);
	     }

           sum += pow(ScatterRspace->AnalyzePower-data[i].anal[j],2)/
	     pow(data[i].Asigma[j],2);
	 }
    }
  //sum *= 10.; //rjc used for fitting n+48Ca
  //-------------------------------------------------------------------

  
 //  sum /=DegreesOfFreedom;
 //  chidif = sum;
   //fit absorption cross section for protons
   double xsec = ScatterRspace->AbsorptionXsection(); //protons
   double const error = 5.;
   cout << sum << " " << xsec << " " << sum +  pow((xsec-543.)/error,2) << endl;
//   return sum +  pow((xsec-543.)/error,2); //fix
 



  // reaction xsections
  
  for(int i=0;i<NXdata;i++)
    {
      double Ecm = Xdata[i].energyCM;
      double Elab = Xdata[i].energyLab;


 
      // integrate wavefunction and find phaseshifts
      ScatterRspace->getSmatrix(Ecm,Elab);

       double xsec = ScatterRspace->AbsorptionXsection(); //protons

       //remove compound elastic part
       if (Ecm < 15.) xsec -= - ScatterRspace->statistical(ScatterRspace->konst,Ecm);

       sumXsec += pow(((xsec-Xdata[i].xsec)/(Xdata[i].sigma)),2); //RJC

    }
  
  //total cross section
  //sumXsec = sumXsec * 10. ;



  if (NTotXdata > 0 && Zp == 0.)
    {
      for(int i=0;i<NTotXdata;i++)
	{

	  double Ecm = TotXdata[i].energyCM;
	  double Elab = TotXdata[i].energyLab;
 

	  // integrate wavefunctions and find S matrices
	  ScatterRspace->getSmatrix(Ecm,Elab);


	  //fit total cross section for neutrons
	  double xsec = ScatterRspace->TotXsection(); //neutrons
	  sumXsec += pow((xsec-TotXdata[i].xsec)/TotXdata[i].sigma,2); //RJC

	}
    }
  

  double total = 0.;
  double chiPoint[4] = {0.,0.,0.,0.};
  if (DegreesOfFreedom >0)
    {
     chiPoint[0] = sum/DegreesOfFreedom*3.;
     total += chiPoint[0];
    }

  

  if (XsecDegreesOfFreedom > 0) 
    {
     chiPoint[2] = sumXsec/XsecDegreesOfFreedom*100.;//*10.;
     total += chiPoint[2];
    }

  if (bprint)
    {
      cout <<"elastic and analyzing power chisq="<< chiPoint[0] << " "  << "reaction xsection="<< chiPoint[2] << endl; 
    }

//double  total = 0.;
  total += BoundChiSquared();
  cout <<" total="<< total <<endl;

  return total;
}

// bound state chisquared//

//********************************

double reaction::BoundChiSquared()
{  


    double Emax =  Efermi ;//- Wstart * fGap_B; // Energy where continuum begins
    double Emin = -200 + Emax;
    std::vector< double > rmesh = BoundRspace->make_rmesh();
    std::vector< double > rmesh_p = BoundRspace->make_rmesh_point();
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
// used to get differen result generated by fitting program and the domgen, so that thi fitting would fit sth else on the data///
// I plot and compare the potentilas and the hamiltonian generated from ./chisq .....(inputfile) x with ./domgenII 
// I wrote down my results , the real nonlocal Potential was different in tow cases
// /////////////////////////////21 Aug 2013 Notes///////////////////////////////////////////////
//testing the hamiltonian
//
//  

/*    std::string hamil_reaction="pot_test/hamil_reaction.out";
    std::ofstream hamilreaction(hamil_reaction.c_str() );  

    for (int ii = 0 ;ii<rmesh.size();++ii){
         for (int kk= 0 ; kk< rmesh.size(); ++kk) { 
              hamilreaction<<rmesh[ii]<<"  "<< rmesh[kk] << " " << BoundRspace->re_hamiltonian(rmesh,-30,0,.5)(ii,kk)<<std::endl;
         }
     }

  //potential test
    std::string fyek_reaction="pot_test/file1_reaction.out";
    std::ofstream fyekreaction(fyek_reaction.c_str() );  

   Pot->setAM(0,.5);
   std::cout<< "V_beta Above "<< Pot->VolumeAbove.beta << " "<<"V_beta_below" <<Pot->VolumeBelow.beta<<std::endl;

   std::cout<<"HF beta0 = "<< Pot->HartreeFock.beta0  << " "<<"HF beta1 = "<< Pot->HartreeFock.beta1  <<std::endl;
   Pot->setEnergy(-30);
   for (int makann = 0 ; makann<rmesh.size();++makann){

   fyekreaction <<rmesh[makann]<<" "<< imag(Pot->nonlocalPart(rmesh[makann],rmesh[makann]))
                               <<" "<< imag( Pot->localPart(rmesh[makann]) )          
                               <<" "<< real( Pot->nonlocalPart(rmesh[makann],rmesh[makann] ))          
                               <<" "<< real( Pot->localPart(rmesh[makann]) )          
   			       <<" "<< Pot->nonlocalIM(rmesh[makann],rmesh[makann])
  			       <<" "<< Pot->nonlocalHF(rmesh[makann],rmesh[makann])
                               <<" "<< Pot->localHF(rmesh[makann])
                               <<" "<< Pot-> HartreeFock.U1(rmesh[makann],rmesh[makann])
                               <<" "<< Pot-> HartreeFock.U0(rmesh[makann],rmesh[makann])
                               <<" "<< real(Pot-> SurfaceAbove.U(rmesh[makann],rmesh[makann]))
                               <<" "<< real(Pot-> SurfaceBelow.U(rmesh[makann],rmesh[makann]))
                               <<" "<< real(Pot-> VolumeAbove.U(rmesh[makann],rmesh[makann]))
                               <<" "<< real(Pot-> VolumeBelow.U(rmesh[makann],rmesh[makann]))
                               <<" "<< imag(Pot-> SurfaceAbove.U(rmesh[makann],rmesh[makann]))
                               <<" "<< imag(Pot-> SurfaceBelow.U(rmesh[makann],rmesh[makann]))
                               <<" "<< imag(Pot-> VolumeAbove.U(rmesh[makann],rmesh[makann]))
                               <<" "<< imag(Pot-> VolumeBelow.U(rmesh[makann],rmesh[makann]))<< std::endl;

}
fyekreaction.close();


   cout<< "S_beta Above " << Pot->SurfaceAbove.betas << endl;
   cout<< "S_beta Belowe "<< Pot->SurfaceBelow.betas << endl;
   cout<< "a Belowe "<< Pot->SurfaceBelow.a << endl;
   cout<< "a Above "<< Pot->SurfaceAbove.a << endl;
   cout<< "R Belowe "<< Pot->SurfaceBelow.R << endl;
   cout<< "R Above "<< Pot->SurfaceAbove.R << endl;

   cout<< "A Belowe "<< Pot->SurfaceBelow.A << endl;
   cout<< "A Above "<< Pot->SurfaceAbove.A << endl;
   cout<< "B Belowe "<< Pot->SurfaceBelow.B << endl;
   cout<< "B Above "<< Pot->SurfaceAbove.B << endl;
   cout<< "D Belowe "<< Pot->SurfaceBelow.C << endl;
   cout<< "D Above "<< Pot->SurfaceAbove.C << endl;

   cout<< "D Belowe "<< Pot->SurfaceBelow.D << endl;
   cout<< "D Above "<< Pot->SurfaceAbove.D << endl;
   cout<< "Wb "<< Pot->SurfaceBelow.Wstart << endl;
   cout<< "W "<< Pot->SurfaceAbove.Wstart << endl;

std::cout<<"DONE the POT test"<<std::endl;

*/

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////End of the test///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



    // Get bound levels
    double tol = 0.01;
    int L_test = 0;
    double J_test = 0.5;
    int N_orbit=0;
    std::vector< lj_eigen_t > bound_levels = 
        BoundRspace->get_bound_levels( rmesh, tol );
    // Create energy meshes
    std::vector< mesh_t > emesh_vec = 
        BoundRspace->get_emeshes( rmesh, Emin, Emax, bound_levels );

    // Calculate propagators
    std::vector< prop_t > prop_vec = 
        BoundRspace->get_propagators( rmesh, emesh_vec );
   
    /// Calculates the Width of 0s1/2//  
    eigen_t bound_info = BoundRspace->find_boundstate( rmesh, Efermi, N_orbit,
                                                       L_test, J_test, tol );

    std::vector< double > e_values = 
        BoundRspace->get_lj_energy_vector( L_test, J_test, emesh_vec );
    
   std::vector< double > s_of_E_QH =  
        BoundRspace->spectral_strength_QH( L_test, J_test, rmesh, emesh_vec,
                                           prop_vec, bound_info.second );

    double gamma1 = BoundRspace->find_width( e_values, s_of_E_QH );
// Calculating the mid point of charge distribution and Particle number as well as R_rms from Point_charge distribution function//
    std::vector<double> PointDist = 
        BoundRspace->point_distribution( rmesh_p, Emax, emesh_vec, prop_vec, 
                                         bound_levels);
//std::cout << "PointDist at r = " << rmesh[0] << "is : "<< PointDist[0]<<" "<<"Emx = "<<Emax << std::endl;
////
//the charge distribution points getting ready to be fit 
/////
////
    double middlehalf=0.476667;// .5;
    double middleone=1.01;
    double middleonehalf=1.47667;//1.5;
    double middle=2.01;
    double middle2half=2.54;
 

    double middlethree=3.0;
    double middlethree1=3.48;
    double middlefour=4.0;	
    double middlefour1=4.48;	
    double middlefive=5.0;	
    double middlefive1=5.48;	
/////
////the unfolded experimental Point(charge) distribution for 14 points(the file: unfolded.out)  :
//for r=01,.47,1.01,1.47,2.01 .... 
//
//
    double unfld[15] = { .0924699 , .0896617 , .0841989 ,
			 .081314 , .0791672, .0739352 , 
			 .0636165 , .0475881 , .0275904 , 
			 .0141982 , .00548428 , .00204393 , 
			 .000574189 , .00019305 };

//Now we should find where the above point are
//since the rmesh in the unfolded exp data is diff from what we are generationg so we have to find the
//closest rmesh[i] to the above rmeshes in the "unfolded.out" then assign the density there to the density we are generationg  and then fit
//well, it is not the best way to do it, it should be optimized, we need to do it once not in each iteration
    int half_index=0;
    int one_index=0;
    int onehalf_index=0;
    int middle_index=0;
    int middle2half_index=0;

    int three_index=0;
    int three1_index=0;
    int four_index=0;
    int four1_index=0;
    int five_index=0;
    int five1_index=0;
//////
    double ParticleNumber=0;
    double minnhalf=100.0;
    double minnone=100.0;
    double minnonehalf=100.0;
    double minn=100.0; 			 
    double minn2half=100.0; 			 
    double minnthree=100.0;
    double minnthree1=100.0;
    double minnfour=100.0; 
    double minnfour1=100.0; 
    double minnfive=100.0; 
    double minnfive1=100.0; 

    std::vector< double > rweights;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        rweights.push_back( BoundRspace->rdelt );
    }

    double pParticleNumber = 0;
    double nParticleNumber = 0; 
    std::vector<double> proton_folding = 
        folded_ch_density( rmesh, rweights, PointDist, 0.5, A );

    std::vector<double> neutron_folding = 
        folded_ch_density( rmesh, rweights, PointDist, -0.5, A );

//finding indexes for rgrids 
//

    for (int i=0; i < rmesh.size() ;++i)
	{      
		if (std::abs(rmesh[i]-middlehalf)<minnhalf) 
			{ half_index=i; 
			  minnhalf=std::abs(rmesh[i]-middlehalf); }

		if (std::abs(rmesh[i]-middleone)<minnone) 
			{ one_index=i; 
			  minnone=std::abs(rmesh[i]-middleone); }
		
		if (std::abs(rmesh[i]-middleonehalf)<minnonehalf) 
			{ onehalf_index=i; 
			  minnonehalf=std::abs(rmesh[i]-middleonehalf); }

		if (std::abs(rmesh[i]-middle)<minn) 
			{ middle_index=i; 
			  minn=std::abs(rmesh[i]-middle); }

		if (std::abs(rmesh[i]-middle2half)<minn2half) 
			{ middle2half_index=i; 
			  minn2half=std::abs(rmesh[i]-middle2half); }

		if (std::abs(rmesh[i]-middlethree)<minnthree) 
			{ three_index=i; 
			  minnthree=std::abs(rmesh[i]-middlethree); }
		
		if (std::abs(rmesh[i]-middlethree1)<minnthree1) 
			{ three1_index=i; 
			  minnthree1=std::abs(rmesh[i]-middlethree1); }

		if (std::abs(rmesh[i]-middlefour)<minnfour) 
			{ four_index=i; 
			  minnfour=std::abs(rmesh[i]-middlefour); }

		if (std::abs(rmesh[i]-middlefour1)<minnfour1) 
			{ four1_index=i; 
			  minnfour1=std::abs(rmesh[i]-middlefour1); }

		if (std::abs(rmesh[i]-middlefive)<minnfive) 
			{ five_index=i; 
			  minnfive=std::abs(rmesh[i]-middlefive); }

		if (std::abs(rmesh[i]-middlefive1)<minnfive1) 
			{ five1_index=i; 
			  minnfive1=std::abs(rmesh[i]-middlefive1); }

        ParticleNumber += 4 * M_PI * PointDist[i] * std::pow(rmesh[i],2) 
                        * rweights[i]; 
        
        pParticleNumber += 4 * M_PI * proton_folding[i] * std::pow(rmesh[i],2) 
                        * rweights[i]; 
        nParticleNumber += 4 * M_PI * neutron_folding[i] * std::pow(rmesh[i],2) 
                        * rweights[i]; 
	}

    std::vector< double > chdf; // folded charge density
    std::vector< double > chdfnew; // folded charge density
    double summs = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double normednew = ( proton_folding[i] + neutron_folding[i] ); 

        double normed = Z * ( proton_folding[i] + neutron_folding[i] ) 
                      / ParticleNumber;

        chdf.push_back( normed );
        chdfnew.push_back( normednew );
		summs += normed * 4 * M_PI * std::pow(rmesh[i],4) * rweights[i];
      }
    double R_rms=std::sqrt( summs / Z ); // Root-mean-square radius
     
    //Calculating Bound ChiSquared
   /* double ch_den_chisq = std::pow( (chdf[0] - 0.0867255)
                        /(.1 * (0.000867255)), 2 );

    double ch_den_half_chisq = std::pow( ( chdf[half_index] - 0.085527 )
                        / (0.000855276), 2 );

    double ch_den_one_chisq = std::pow( ( chdf[one_index] - 0.0825104)
                        / ( 0.000825104 ), 2 );

    double ch_den_onehalf_chisq = std::pow( ( chdf[onehalf_index] - 0.0800167)
                        / (0.000800167), 2 );

    double ch_den_middle_chisq = std::pow( ( chdf[middle_index] - 0.0770948 ) 
                        /(.000770948), 2 );

    double ch_den_middle2half_chisq = std::pow( ( chdf[middle2half_index] - 0.0709657) 
                        /(0.000709657), 2 );

    double ch_den_three_chisq = std::pow( ( chdf[three_index] - 0.0605082 )
                        / (  0.000605082), 2 );

    double ch_den_three1_chisq = std::pow( ( chdf[three1_index] - 0.0457602 )
                        / (  0.000457602), 2 );

    double ch_den_four_chisq = 2. * std::pow( ( chdf[four_index] - 0.0278664 )
                        / (  0.000278664), 2 );
    
    
    double ch_den_four1_chisq = .5 * std::pow( ( chdf[four1_index] - 0.0153698 )
                        / (  0.000153698), 2 );
    
    double ch_den_five_chisq = .02 * std::pow( ( chdf[five_index] - 0.00658348 )
                        / (  0.0000658348), 2 );
    

    double ch_den_five1_chisq = .005 * std::pow( ( chdf[five1_index] - 0.00274881 )
                        / (  0.0000274881), 2 );
*/

//
//the UNFOLDED charge distribution is being fitted here 
//the unfolded experimental data are in the unfolded.out
//
    double ch_den_chisq = std::pow( (PointDist[0] -unfld[0]) 
                        /( (unfld[0]/1000.0)), 2 );

    double ch_den_half_chisq = std::pow( ( PointDist[half_index] - unfld[1] )
                        / (unfld[1]/1000.0), 2 );

    double ch_den_one_chisq = std::pow( ( PointDist[one_index] - unfld[2])
                        / (unfld[2]/100.0 ), 2 );

    double ch_den_onehalf_chisq = std::pow( ( PointDist[onehalf_index] - unfld[3])
                        / (unfld[3]/100.0), 2 );

    double ch_den_middle_chisq = std::pow( ( PointDist[middle_index] - unfld[4] ) 
                        /(unfld[4]/100.0), 2 );

    double ch_den_middle2half_chisq = std::pow( ( PointDist[middle2half_index] - unfld[5]) 
                        /(unfld[5]/100.), 2 );

    double ch_den_three_chisq = std::pow( ( PointDist[three_index] - unfld[6] )
                        / ( unfld[6]/100.0), 2 );

    double ch_den_three1_chisq = std::pow( ( PointDist[three1_index] - unfld[7] )
                        / ( unfld[7]/100.0), 2 );

    double ch_den_four_chisq = 2. * std::pow( ( PointDist[four_index] - unfld[8] )
                        / ( unfld[8]/100.0), 2 );
    
    
    double ch_den_four1_chisq = .5 * std::pow( ( PointDist[four1_index] - unfld[9] )
                        / (  unfld[9]/100.0), 2 );
    
    double ch_den_five_chisq = .02 * std::pow( ( PointDist[five_index] - unfld[10] )
                        / (  unfld[10]/100.0), 2 );
    

    double ch_den_five1_chisq = .005 * std::pow( ( PointDist[five1_index] - unfld[11] )
                        / (  unfld[11]/100.0), 2 );

    double R_rms_chisq = std::pow( ( R_rms - Chi_Squared_data.R_rms_exp )
                       /(.2 * Chi_Squared_data.err_R_rms_exp), 2 );

    double gamma1_chisq = std::pow( ( gamma1 - 21.3 ) / 0.9, 2 );	
  
    
   /// Energy levels Chisq
   cout<<Nlevel<<endl;
   double Energy_chisq = 0;
   for (int i = 0 ; i < Nlevel ; ++i) { 
	
		
        if (Level[i].Efit == 1) {
		eigen_t Th_Energy = BoundRspace->find_boundstate( rmesh, Efermi, Level[i].N,
                                            Level[i].l ,  Level[i].j , tol );

                std::cout<<i<<" "<<"theory E is = "<< " " << Th_Energy.first<< "\t "<< " exp Energy is = " << " "<< Level[i].Energy<<std::endl; 
		Energy_chisq += std::pow( ( Th_Energy.first - Level[i].Energy ) 
                      / (3.0 * Level[i].SigmaEnergy ), 2 );
	    }
   }	
   
    ///////

    // determine experimental particle number (make it a little less 
    // to decrease the chance of the calculated number going over the
    // experimental one)
    double exp_particle_number;
    double err_particle_number = 0.02;
    if ( Zp == 1 ) exp_particle_number = Z - 0.01;
    else exp_particle_number = A - Z - 0.01;

    double particle_chisq = std::pow((( ParticleNumber - exp_particle_number ) 
                          /(3.  * err_particle_number)), 2 );
   
    // Add all the pieces together
    //
   double BoundChisquared =  tot_chisq_p_250 + tot_chisq_p_650;

//   double BoundChisquared =  ch_den_chisq +  ch_den_half_chisq + ch_den_one_chisq + 
 //         ch_den_onehalf_chisq + ch_den_middle_chisq + ch_den_middle2half_chisq + Energy_chisq; 
 //  double BoundChisquared =  ch_den_chisq + ch_den_five1_chisq+R_rms_chisq + particle_chisq ; 
//   double BoundChisquared =  ch_den_five1_chisq ; 
   

//////
//Printing out the results to check and observe
//////
//
    cout << "charge density at 0 =  " << PointDist[0] << " " 
         << "chisq = " << ch_den_chisq << " "<< "exp = "<< unfld[0] << endl;

    cout << "charge density at 0.5 =  " << PointDist[half_index] << " " 
         << "chisq = "  << ch_den_half_chisq << " "<< "exp = "<< unfld[1] << endl;

    cout << "charge density at 1 =  " << PointDist[one_index] << " " 
         << "chisq = "  << ch_den_one_chisq << " "<< "exp = "<< unfld[2] << endl;

    cout << "charge density at 1.5 =  " << PointDist[onehalf_index] << " " 
         << "chisq = "  << ch_den_onehalf_chisq << " "<< "exp = "<< unfld[3] << endl;

    cout << "charge density at 2 =  " << PointDist[middle_index] << " " 
         << "chisq = "  << ch_den_middle_chisq << " "<< "exp = "<< unfld[4]  << endl;

    cout << "charge density at 2.5 =  " << PointDist[middle2half_index] << " " 
         << "chisq = "  << ch_den_middle2half_chisq << " "<< "exp = "<< unfld[5]  << endl;

    cout << " charge density at 3 = " << PointDist[three_index] << " " 
         << "chisq  = "  << ch_den_three_chisq << " "<< "exp = "<< unfld[6] <<endl;

    cout << " charge density at 3.5 = " << PointDist[three1_index] << " " 
         << "chisq  = "  << ch_den_three1_chisq << " "<< "exp = "<< unfld[7]  <<endl;

    cout << " charge density at 4 = " << PointDist[four_index] << " " 
         << "chisq  = " << ch_den_four_chisq  << " "<< "exp = "<< unfld[8] <<endl;

    cout << " charge density at 4.5 = " << PointDist[four1_index] << " " 
         << "chisq  = " << ch_den_four1_chisq  << " "<< "exp = "<< unfld[9]  <<endl;

    cout << " charge density at 5 = " << PointDist[five_index] << " " 
         << "chisq  = " << ch_den_five_chisq  << " "<< "exp = "<< unfld[10] <<endl;

    cout << " charge density at 5.5 = " << PointDist[five1_index] << " " 
         << "chisq  = " << ch_den_five1_chisq  << " "<< "exp = "<< unfld[11]  <<endl;

    cout << Nlevel <<  " " << "Energy chisq = " << Energy_chisq << endl;
    cout << " " <<endl;    
    cout << " " <<endl;    
    cout << "R_rms = " << R_rms  << " " 
         << "R_rms_chisq = " << R_rms_chisq << endl ;

    cout << " " <<endl;    
    cout << "particle number = " << ParticleNumber << " " 
         << "particle_chisq = " << particle_chisq << endl ;

    cout << " " <<endl;    
    cout << "gamma1 = "<< gamma1 << " " 
         << "gamma1_chisq = " << gamma1_chisq << endl ;

    cout << " " <<endl;    
    cout << "total bound chisq = " << BoundChisquared << endl ;
    cout << " " <<endl;    
   //p(n)ParticleNumber is from folding 
    

/////
//////saving the Point Dist in a file, it could be useful for checking especially in the first iteration, I used it
/////
//
    std::string akbar="akbarr.out";
    std::ofstream filedist(akbar.c_str() );  
    for (int kk=0; kk < rmesh.size() ; kk++) {
         filedist << rmesh[kk] << " " << PointDist[kk] <<std::endl;
     }
    filedist.close();	 
 //   cout << PointDist[0] << " " << proton_folding[0] << " " << neutron_folding[0] << " " <<  nParticleNumber << " " << pParticleNumber << endl ;
    return BoundChisquared;
}


//*********************************************************************
  /**
   *relativistic conversion from Lab kinetic energyto total CM kinetic energy
   \param Elab energ in laboratory frame [MeV]
  */
double reaction::energyLab2Cm(double Elab)
{
  if (Elab < 0.) return Elab*A/(1.+A);
 // center of mass velocity in units of c
 double vcm = sqrt(Elab*(Elab+2.*scatterRspace::m0))/(Elab+(1.+A)*
 scatterRspace::m0);
 //gamma factor for this velocity
 double gam = 1./sqrt(1.-std::pow(vcm,2));
 double Ecm = (gam-1.)*(1.+A)*scatterRspace::m0 +
               gam*Elab*(A-1.)*scatterRspace::m0/((A+1.)*
              scatterRspace::m0+Elab);
 return Ecm;
}
//**********************************************************************
  /**
   *relativistic conversion from Lab kinetic energy to total CM kinetic energy
   \param Ecm is energy in the center-of-mass frame
  */
double reaction::energyCm2Lab(double Ecm)
{
  if (Ecm < 0.) return Ecm/A*(1.+A);

  //find momentum of projecile, also momentum of target
  double pc = sqrt(Ecm)*sqrt((Ecm+2.*scatterRspace::m0)*(Ecm+2.*A*
              scatterRspace::m0)*
              (Ecm+2.*(A+1.)*scatterRspace::m0))/2./(Ecm+(A+1)
               *scatterRspace::m0);
  //velocity of target in units of c
  double vtarget = pc/sqrt(pow(pc,2)+pow(A*scatterRspace::m0,2));
  //gamma factor for this velocity
  double gam = 1./sqrt(1.-pow(vtarget,2));
  // tot energy of projectile (neutron or proton in com frame)
  double Eproj = sqrt(pow(scatterRspace::m0,2)+pow(pc,2));
  double Elab = (Eproj + vtarget*pc)*gam;
  //this energy contains rest mass , so remove it 
  Elab -= scatterRspace::m0;
  return Elab;
}

void reaction::Print_DiffXsection( double Elab ) {

    double Ecm = energyLab2Cm( Elab );
    ScatterRspace->getSmatrix(Ecm,Elab);

    double theta_min = 10;
    double theta_max = 180;
    double theta_pts = 100;
    double delta_theta = ( theta_max - theta_min ) / theta_pts;

    ofstream file( "diff_cross_section.out" );
    for ( int i = 0; i < theta_pts; i++) {

        double theta_deg = theta_min + i * delta_theta; // angle in degrees
	    double theta_rad = theta_deg * pi/180.; // angle in radians
        file << theta_deg << " " 
             <<  ScatterRspace->DifferentialXsection( theta_rad ) << endl;

    }

    file.close();
    file.clear();
}





