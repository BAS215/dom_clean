#include "ratio.h"
/**
 * Constructor
\param title0 gives the root filename for input data
\param Reaction10 specifies reaction 1
\param Reaction20 specifies reaction 2
*/

ratio::ratio(string title0,reaction* Reaction10, reaction* Reaction20)
{

  title = title0;

  Reaction1 = Reaction10;
  Reaction2 = Reaction20;

  string directory("");

  //open data file
  string filename = directory + title + ".data";
  cout << filename << endl;
  ifstream fileData (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileData.fail()) 
    {
      cout << "couldn't open data file" << fileData << endl;
      abort();
    }

  else
    {
      string line;
      getline(fileData,line);
      int NN;
      fileData >> NN;
      Rdata = new rdata [NN];
      Ndata = 0;
      for (int i=0;i<NN;i++)
	{
	  fileData >> Rdata[Ndata].energyLab >> Rdata[Ndata].x >> 
          Rdata[Ndata].sigma;
          if (Rdata[Ndata].energyLab < 5.) continue;
          //Rdata[Ndata].x /= 200.;
          //Rdata[Ndata].sigma /= 200.;

          Rdata[Ndata].energyCM1 = Reaction1->
                              energyLab2Cm(Rdata[Ndata].energyLab);
          Rdata[Ndata].energyCM2 = Reaction2->
                              energyLab2Cm(Rdata[Ndata].energyLab);
          Ndata++;
	}
      fileData.close();
      fileData.clear();
    }
  cout << "Ndata = " << Ndata << endl;
}
//******************************************************
  /**
   *Destructor
   */
ratio::~ratio()
{
  delete [] Rdata;
}
//********************************************************
  /**
   * calculates the ratio of total cross sections, i.e.,
   *\f$ \frac{\sigma_2-sigma_1}{\sigma_2+\sigma_1} \f$
   \param Elab is the laboratory energy in MeV
   */ 
double ratio::getRatio(double Elab)
{
  double Ecm = Reaction1->energyLab2Cm(Elab);
  Reaction1->InitializeForEcm(Ecm,Elab);
  Ecm = Reaction2->energyLab2Cm(Elab);
  Reaction2->InitializeForEcm(Ecm,Elab);

  //integrate wave functions and find phase shift for both reactions
  Reaction1->scatter->integrateWave();
  Reaction2->scatter->integrateWave();

  double x1 = Reaction1->scatter->TotXsection();
  double x2 = Reaction2->scatter->TotXsection();
  return (x2-x1)/(x2+x1); 

}

//******************************************************
  /**
   * returns the chi squared per degree of freedom between the 
   * calculated ratios and the experimental data
   */
double ratio::ChiSquared()
{



  double sum = 0.;
  for (int i=0;i<Ndata;i++)
    {
      double Elab = Rdata[i].energyLab;
      double Ecm = Rdata[i].energyCM1;
      Reaction1->InitializeForEcm(Ecm,Elab);
      Ecm = Rdata[i].energyCM2;
      Reaction2->InitializeForEcm(Ecm,Elab);
  

      //integrate wave functions and find phase shift for both reactions
      Reaction1->scatter->integrateWave();
      Reaction2->scatter->integrateWave();


      double x1 = Reaction1->scatter->TotXsection();
      double x2 = Reaction2->scatter->TotXsection();
      double x = (x2-x1)/(x2+x1); 
      sum += pow((x-Rdata[i].x)/Rdata[i].sigma,2);
   }

  sum /= (double)Ndata;

  sum *= 8.;
  cout << "ratio " << sum << endl;


  return sum;
}
//***********************************************************
  /**
   * Make a root plot of the experimental ratio and the 
   * calculated values
   */
void ratio::PlotFit()
{
#ifdef root
  TCanvas *canvas = new TCanvas("fit");
  canvas->SetLogx();
  TH2S * frame = new TH2S("frame","",10,5,300,10,-.02,.15);
  frame->GetXaxis()->SetTitle("E_{lab}");
  frame->GetYaxis()->
         SetTitle("(#sigma_{2}-#sigma_{1})/(#sigma_{2}+#sigma_{1})");
  frame->SetStats(kFALSE);
  frame->Draw();

  ofstream xFile("ratiox.dat");
  xFile << Ndata << endl;
  if (Ndata > 0)
    {
     double xExp[Ndata];
     double yExp[Ndata];
     double xErrorExp[Ndata];
     double yErrorExp[Ndata];

     for (int i=0;i<Ndata;i++)
       {
         xExp[i] = Rdata[i].energyLab;
         yExp[i] = Rdata[i].x;
         xErrorExp[i] = 0.;
         yErrorExp[i] = Rdata[i].sigma;
         xFile << xExp[i] << " " << yExp[i] << " " << yErrorExp[i] << endl;
       }

     TGraphErrors * graphExp = new TGraphErrors(Ndata,xExp,yExp,xErrorExp,
                                  yErrorExp);
     graphExp->SetMarkerStyle(21);
     graphExp->Draw("P");
    }
  xFile.close();
  ofstream ofFile("ratio.dat");
  int const Npoints = 100;
  double xFit[Npoints];
  double yFit[Npoints];
  double xStart = log(5.);
  double xStop = log(300);
  for (int i=0;i<Npoints;i++)
    {
      double xx = (xStop-xStart)/((double)Npoints)*(double)i + xStart;
      double Elab = exp(xx);
      xFit[i] = Elab;
      yFit[i] = getRatio(Elab);
      ofFile << xFit[i] << " " << yFit[i] << endl;
    }
  TGraph * graph = new TGraph(Npoints,xFit,yFit);
  graph->Draw("L");

  canvas->Write();

  ofFile.close();
  ofFile.clear();


  

#endif
}
//**************************************************************************
  /**
   * opens root file to save spectra in
   */
void ratio::OpenRootFile()
{
#ifdef root
  string filename(title + ".root");
  cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);

  f = new TFile(filename.c_str(),"RECREATE");
#endif
}
//**************************************************************************
  /**
   * closes root file
   */
void ratio::CloseRootFile()
{
#ifdef root
  f->Write();
  f->Close();
#endif
}

