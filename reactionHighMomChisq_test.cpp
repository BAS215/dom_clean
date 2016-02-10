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
  Pot->load( Rc, VHFvol, VHFsur, RHF, aHF,  RHFs, aHFs, beta_nl_R0, AHF, beta_nl_R1, RsurfaceAbove, RsurfaceBelow, asurfaceAbove, asurfaceBelow, 
	     Asurface_A, Asurface_B, BsurfaceA, CsurfaceA, DsurfaceA, Bsurface, Csurface, Dsurface,
	      Wstart*fGap_A, Wstart*fGap_B, Efermi, beta_nl_I0, beta_nl_I1, beta_nl_I0_sur, beta_nl_I1_sur,
             RvolumeAbove,RvolumeBelow, deltaRvolume, expRvolume, avolumeAbove,avolumeBelow, Avolume_A, Avolume_B, BvolumeAbove, BvolumeBelow , EpvolumeAbove,EpvolumeBelow,
	     mvolume, AsyVolume, alphaVolume, EaVolume_a, EaVolume_b, Rso, aso, Vso, AWso, BWso ,
	     V_wine , R_wine , rho_wine);

  return true;
  //return AdjustFermi();
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


     // cout <<"elastic and analyzing power chisq="<< chiPoint[0] << " " << chiPoint[1] << " " << "reaction xsection="<< chiPoint[2] << endl; 
  total = 0.;
  total += BoundChiSquared();
  cout <<" total="<< total <<endl;

  return total;
}

// bound state chisquared//

//********************************

double reaction::BoundChiSquared()
{  


    double Emax =  Efermi - Wstart * fGap_B; // Energy where continuum begins
    double Emin = -200 + Emax;
    std::vector< double > rmesh = BoundRspace->make_rmesh();

    // Get bound levels
    double tol = 0.05;
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
    ////
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//********************************alculating the missing S(E,P) Chi Squared***************************************************************
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 


    std::string input_dir = "Input/";
    int num_lj = 11;
    std::string p_filename = input_dir + "pca40.inp";

    // Create Nuclear Parameter Objects
    NuclearParameters Nu_p = read_nucleus_parameters( p_filename );

    std::vector< NuclearParameters > Nu_vec;
    Nu_vec.push_back( Nu_p );

    std::string parameters_filename = "hosfit.inp"; 
 
    // Read in DOM parameters
    std::ifstream pfile( parameters_filename.c_str() );
    if ( pfile.is_open() !=1 ) {
       	  std::cout << "could not open file " << parameters_filename << std::endl;
       	  std::abort();
    }
    pfile.close();
    pfile.clear();
      
    int lmax = 5;
    int nu=0;
    double tz = .5 ; //nu - 0.5; // +0.5 for protons, -0.5 for neutrons
    double calc_norm = 20 ; 
    // Create momentum space grid for E-slice
    std::vector<double> kmesh2; // wave number in units of fm^{-1}
    std::vector<double> pmesh; // momentum in units of MeV / c
    double pmin = 250; 
    double pmax = 650;
    double deltap = 80;
    int ppts = static_cast<int> ( ( pmax - pmin ) / deltap ) + 1;
    double p=0;
    for ( int i = 0; i < ppts; ++i ) {

         	p = pmin + i * deltap;
       	 	pmesh.push_back( p );

       	         double k = p / hbarc;
                kmesh2.push_back( k );
    }

    // Nucleus Object
    const NuclearParameters &Nu = Nu_vec[0];

    // Construct Parameters Object
    Parameters pap = get_parameters( parameters_filename, Nu.A, Nu.Z,1 );

    // Construct Potential Object
    // type=1
    pot U = get_bobs_pot2( 1 , 4, 1, tz, Nu, pap );
    pot *U1 = &U;


    // Construct Object for calculating bound-state properties
    boundRspace B(12.0, 180 , Nu.Ef, Nu.ph_gap, lmax, Nu.Z, 1, Nu.A, U1 );
    // Create Energy grid for E-slice
    double Emin2 = -300;
    double Emax2 = -25;
    double deltaE2 = 2;
    int epts2 = static_cast<int>( ( Emax2 - Emin2 ) / deltaE2 ) + 1;
    std::vector<double> emesh2;
    for ( int i = 0; i < epts2; ++i ) {

            emesh2.push_back( Emin2 + i * deltaE2 );
    }

    matrix_t S_of_kE_mtx2( kmesh2.size(), emesh2.size() );

    // intialize matrices to zero

    for ( unsigned int i = 0; i < kmesh2.size(); ++i ) {
         for ( unsigned int j = 0; j < emesh2.size(); ++j ) {

              S_of_kE_mtx2( i, j ) = 0;
         }
    }

   // Loop over lj channels
   for ( int L = 0; L < lmax + 1; ++L ) {

        for( int up = -1; up < 2; up+=2 ) {

                double xj = L + up / 2.0;
                int j_index = ( up - 1 ) / 2;
                if ( xj < 0 ) continue;



                // Create Bessel Function matrix in k and r
                matrix_t bess_mtx2( kmesh2.size(), rmesh.size() );
                for( unsigned int nk = 0; nk < kmesh2.size(); ++nk ) {
                for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

                    double rho = kmesh2[nk] * rmesh[nr];

                    bess_mtx2( nk, nr ) = gsl_sf_bessel_jl( L, rho );
                }
                }

                // Calculate S( k; E ) for E-slice 

                for ( unsigned int m = 0; m < emesh2.size(); ++m ) {

            
                    double E = emesh2[m];

                    // Propagator
                    cmatrix_t G = BoundRspace->propagator( rmesh, E, L, xj );
            //        if (L==0 && m==0){std::cout<<"G(0,0)="<< G(0,0)<<std::endl;}  
                    // Spectral Function in momentum space, S( k; E ) 
                    for( unsigned int nk = 0; nk < kmesh2.size(); ++nk ) {

                        double rsums = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                            double jl1 = bess_mtx2( nk, i );

                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            double jl2 = bess_mtx2( nk, j );
                    
                            rsums -= rmesh[i] * jl1 * imag( G( i, j ) ) 
                                   * rmesh[j] * jl2 *BoundRspace->rdelt * 2 / M_PI / M_PI;
                        }
                        } // end loop over radial coordinates
                        S_of_kE_mtx2( nk, m ) += ( 2 * xj + 1 ) * rsums;
				
                    } // end loop over k

               } // end loop over energy


           }// end loop over j

     } // end loop over L
       
     // write in units of MeV^-4 sr^-1
     double fac = std::pow( hbarc, 3 ) * 4 * M_PI;

//    S(E_m,P_m) Chi squared COnstruction

     
     int point250 [] = {124 , 118 , 113, 106 , 101 , 95 , 89 , 83 , 77};
 
     int point650 [] = {119 , 106 , 93, 81, 67 , 53 , 41};

     double point250_data [] = {1.6803e-11 , 1.2703e-11, 9.7322e-12, 7.212e-12 ,
				   5.5253e-12 , 4.2899e-12, 3.3753e-12, 2.7274e-12, 2.3556e-12 };

     double point650_data [] = {1.6746e-14 , 2.581e-14 , 2.7952e-14 , 2.5291e-14 , 2.7942e-14 , 2.2278e-14 , 2.3966e-14  };

     double err_250 = 1e-12; 
     double err_650 = 0.5e-12; 

     std::vector <double> chisq_p_250;
     std::vector <double> S_missing_250;
     double tot_chisq_p_250 = 0; 

     std::vector <double> chisq_p_650;
     std::vector <double> S_missing_650;
     double tot_chisq_p_650 = 0; 

     for  (int i = 0 ; i < 9 ; ++i) {
                S_missing_250.push_back(S_of_kE_mtx2( 0, point250[i]-1 ) / fac / calc_norm);

		chisq_p_250.push_back(pow(((S_missing_250[i]-point250_data[i]) / err_250) , 2));

		tot_chisq_p_250 += chisq_p_250[i];

//             std::cout << chisq_p_250[i] << " " << "calculated 250"<< " " << i << " = "
//		          << S_missing_250[i] << " " << "Daniela 250 " << i << " = " 
//	                  << point250_data[i] <<std::endl;
     }

     for  (int i = 0 ; i < 7 ; ++i) {
                S_missing_650.push_back(S_of_kE_mtx2( 0, point650[i]-1 ) / fac / calc_norm);

		chisq_p_650.push_back(pow(((S_missing_650[i]-point650_data[i]) / err_650) , 2));

		tot_chisq_p_650 += chisq_p_650[i];

//             std::cout << chisq_p_650[i] << " " << "calculated 650"<< " " << i << " = "
//		          << S_missing_650[i] << " " << "Daniela 650 " << i << " = " 
//	                  << point650_data[i] <<std::endl;
     }
     tot_chisq_p_250 = tot_chisq_p_250 ;	
     std::cout<<"total 250 chisq =  "<<" "<<tot_chisq_p_250<<std::endl;
     tot_chisq_p_650 = tot_chisq_p_650 ;	
     std::cout<<"total 650 chisq =  "<<" "<<tot_chisq_p_650<<std::endl;

//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 
//***************************************************************************************************************************************** 


// Calculating the mid point of charge distribution and Particle number as well as R_rms from Point_charge distribution function//

    // Add all the pieces together
    
    double BoundChisquared = tot_chisq_p_650+tot_chisq_p_250 ;

//    filedist.close();	 
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





