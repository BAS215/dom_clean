#include "reaction.h"

//***********************************************************************
  /**
   *plots the fit to the elastic scattering angular distributions
   * and analyzing powers- output is wriiten to root file title.root
   */
void reaction::PlotFit()
{

#ifdef root
  double factA = 0;

  TCanvas *fit = new TCanvas("fit");
  fit->Divide(3,1);
  TPad* fit_1 = (TPad*)(fit->GetPrimitive("fit_1"));
  fit_1->SetLogy();
  TPad* fit_2 = (TPad*)(fit->GetPrimitive("fit_2"));
  fit_2->SetLogy();

  fit->cd(1);

  // axies
  int ymin = (int)floor(log10(xsecMin)) - Ndata;
  int ymax = (int)ceil(log10(xsecMax));
  TH2S *hist = new TH2S("hist",title.c_str(),10,0,180.
                        ,10,pow(10.,ymin),pow(10.,ymax));
  hist->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
  hist->SetStats(kFALSE);
  hist->Draw();


  TH2S *hist2 = new TH2S("hist2",title.c_str(),10,0,180.,10,.1,pow(8.,Ndata));
  hist2->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist2->GetYaxis()->SetTitle("#sigma/#sigma_{Ruth}");
  hist2->SetStats(kFALSE);
  fit->cd(2);
  hist2->Draw();


  TH2S *hist3 = new TH2S("hist3",title.c_str(),10,0,180.,10,-1.,32.);
  hist3->GetXaxis()->SetTitle("#theta_{c.m.} [deg]");
  hist3->GetYaxis()->SetTitle("Analyze Power");
  hist3->SetStats(kFALSE);
  fit->cd(3);
  hist3->Draw();

  fit->cd(1);

  TGraphErrors *gdata[Ndata];
  TGraph *GA[Ndata];
  TGraph *GAT[Ndata];
  TLine  *LA[Ndata];
  double zero[100];
  for (int i=0;i<100;i++) zero[i] = 0.;

  string filename = title+"x.dat";
  ofstream fdata(filename.c_str());
  filename = title+"xx.dat";
  if (Zp == 1) ofstream fdatax(filename.c_str());

  int jj = 0;
  for (int i=0;i<Ndata;i++)
    {

      if (data[i].nX == 0) continue;
      fdata << data[i].energyLab << " " << data[i].nX <<  endl;

      double y[data[i].nX];
      double error[data[i].nX];
      double fact = pow(10.,jj);
      for (int j=0;j<data[i].nX;j++)
	{
	  y[j] = data[i].xsec[j]/fact;
          error[j] = data[i].Xsigma[j]/fact;
	  fdata << data[i].Xtheta[j] << " " << data[i].xsec[j] << " " <<
	    data[i].Xsigma[j] << endl;
	}
      gdata[i] = new TGraphErrors(data[i].nX,
                data[i].Xtheta,y,zero,error);
      gdata[i]->SetLineColor(i%2+1);
      gdata[i]->SetMarkerStyle(21);
      gdata[i]->SetMarkerSize(.2);
      gdata[i]->SetMarkerColor(i%2+1);
      gdata[i]->Draw("P");
      jj++;
    }
  fdata.close();
  fdata.clear();
  filename = title+"c.dat";
  fdata.open(filename.c_str());
  filename = title+"a.dat";
  ofstream hdata(filename.c_str());


  filename = title+"xx.dat";
  ofstream fxdata(filename.c_str());

  TGraph *Gdata[Ndata];
  TGraph *Gratio[Ndata];
  TGraphErrors *GratioExp[Ndata];
  TLine  *Lratio[Ndata];
  jj = 0;


  prepare();

  for(int i=0;i<Ndata;i++)
    {

      if (data[i].nX > 0)fdata << data[i].energyLab << " " << 177 << endl;
      double Ecm =data[i].energyCM;
      double Elab =data[i].energyLab;




       // integrate wavefunction and find phaseshifts
       ScatterRspace->getSmatrix(Ecm,Elab);
       if (Ecm < 15.)ScatterRspace->statistical(ScatterRspace->konst,Ecm);
       if (data[i].nX > 0)
	 {
          double ang[177];
          double yy[177];
          double rr[177];
          double fact = pow(10.,jj);
          double factr = pow(10.,jj);
          // loop over angles
	  fxdata << 177 << endl;
          for (int j=3;j<180;j++)
	    {
	      double angle = (double)j*pi/180.;
              yy[j-3] = ScatterRspace->DifferentialXsection(angle);
	      if (Ecm < 15.)
		{
                  double extra = ScatterRspace->DifferentialXsectionCE(angle);
		  yy[j-3] += extra;
		}
              ang[j-3] = (double)j;
	      double ratioR = yy[j-3]/ScatterRspace->Rutherford(angle);
              rr[j-3] = ratioR*factr;
              fdata << ang[j-3] << " " << yy[j-3] << endl;
              yy[j-3] /= fact;

              fxdata << angle << " " << rr[j-3]<< endl;
	    }
          double rrexp[data[i].nX];
          double rrexpSigma[data[i].nX];
          fxdata << data[i].nX << " " << Elab << endl;
          for (int j=0;j<data[i].nX;j++)
	    {
	      double angle = data[i].Xtheta[j]*pi/180.;
              double ruth = ScatterRspace->Rutherford(angle);
              rrexp[j] = data[i].xsec[j]/ruth*factr;
              rrexpSigma[j] = data[i].Xsigma[j]/ruth*factr;
              fxdata << angle << " " << rrexp[j] 
              << " " << rrexpSigma[j] << endl;
   	    }
          Lratio[i] = new TLine(0.,factr,180.,factr);
          Lratio[i] -> SetLineColor(i%2+1);

          Gdata[i] = new TGraph(177,ang,yy);
          Gdata[i]->SetLineColor(i%2+1);
          Gdata[i]->Draw("C");
          Gratio[i] = new TGraph(177,ang,rr);
          Gratio[i]->SetLineColor(i%2+1);
          GratioExp[i] = new TGraphErrors(data[i].nX,data[i].Xtheta,rrexp,
                              zero,rrexpSigma);
          GratioExp[i]->SetLineColor(i%2+1);
          GratioExp[i]->SetMarkerStyle(21);
          GratioExp[i]->SetMarkerSize(.5);
          GratioExp[i]->SetMarkerColor(i%2+1);
          fit->cd(2);
          Lratio[i]->Draw();
          GratioExp[i]->Draw("P");
          Gratio[i]->Draw("C");
          fit->cd(1);
	  jj++;
	 }

       //see if there is analyzing power data
       if (data[i].nA > 0)
	 {
           hdata << data[i].energyLab << " " << data[i].nA << endl;

           double aa[data[i].nA];
	   for (int j=0;j<data[i].nA;j++)
	     {
	       aa[j] = data[i].anal[j] + factA;
	       hdata << data[i].Atheta[j] << " " 
		     << data[i].anal[j] << " " << data[i].Asigma[j] << endl;
	     }
           LA[i] = new TLine(0.,factA,180.,factA);
           LA[i]->SetLineColor(i%2+1);
	   GA[i] = new TGraphErrors(data[i].nA,data[i].Atheta,aa,zero,
           data[i].Asigma);
           GA[i]->SetMarkerStyle(21);
           GA[i]->SetMarkerSize(.5);
           GA[i]->SetMarkerColor(i%2+1);
	   fit->cd(3);
           GA[i]->Draw("P");
           LA[i]->Draw();
	   double aaa[177];
	   double spinR[177];
           double ang[177];
           hdata << 177 << endl;
	   for (int j=3;j<180;j++)
	     {
	       double angle = (double)j*pi/180.;
	       double shape = ScatterRspace->DifferentialXsection(angle);


	       ang[j-3] = (double)j;
	       aaa[j-3] = ScatterRspace->AnalyzePower ;
           
               spinR[j-3] = ScatterRspace->SpinRotation;
               if (Ecm < 15.)
		 {
                  double compound = 
                  ScatterRspace->DifferentialXsectionCE(angle);
                  aaa[j-3] *= shape/(shape+compound);
                  spinR[j-3] *= shape/(shape+compound);
		 }
	       hdata << ang[j-3] << " " << aaa[j-3] << " " 
                     << spinR[j-3] << endl;
               aaa[j-3] += factA;
	     }
           GAT[i] = new TGraph(177,ang,aaa);
           GAT[i]->SetLineColor(i%2+1);
           GAT[i]->Draw("C");
           factA += 2.;
	   fit->cd(1);
	 }
    }
   fdata.close();
   fdata.clear();
   hdata.close();
   hdata.clear();
   fxdata.close();
   fxdata.clear();
   filename = title+"d.dat";
   fdata.open(filename.c_str());

  //plot compound elastic for lowest energy
   double Ecm =data[0].energyCM;
   double Elab =data[0].energyLab;

   fdata << Elab << " " << 177 << endl;


   // integrate wavefunction and find phaseshifts
   ScatterRspace->getSmatrix(Ecm,Elab);
   //get cross section
   double xabs = ScatterRspace->AbsorptionXsection();

   if (Ecm < 30.)
     {
     //calulate CN decay widths
     ScatterRspace->statistical(ScatterRspace->konst,Ecm);

     cout << "absorption + compound elastic = " << xabs << " mb" << endl;
     cout << "absorption = " << ScatterRspace->sigmaAbsorption 
        << " mb " << endl;
     cout << "compound elastic = " << ScatterRspace->sigmaCompoundElastic 
        << " mb " << endl;

     // loop over angles
     double yy[177];
     double zz[177];
     double ang[177];
     for (int j=3;j<180;j++)
        {
         double angle = (double)j*pi/180.;
         yy[j-3] = ScatterRspace->DifferentialXsectionCE(angle);
         zz[j-3] = ScatterRspace->DifferentialXsection(angle)+ yy[j-7];

         ang[j-3] = (double)j;
         double  aaa = ScatterRspace->AnalyzePower ;
         //cout << aaa << endl;
         aaa *= (zz[j-3]-yy[j-3])/zz[j-3];
         fdata << ang[j-3] << " " << yy[j-3] << endl;
         }
      TGraph CNelastic (177,ang,yy);
      CNelastic.SetLineColor(3);
      CNelastic.SetLineStyle(2);
      CNelastic.Draw("C");
      //total
      TGraph totElastic(177,ang,zz);
      totElastic.SetLineStyle(1);
      totElastic.SetLineColor(2);
      totElastic.Draw("C");
     }
   else
     {
       ScatterRspace->sigmaAbsorption = 0.;
       ScatterRspace->sigmaCompoundElastic = 0.;
     }

  fit->Write();
  fdata.close();
  fdata.clear();

#endif
}
//****************************************************************************
  /**
   *plots the fitted potential and many other graphs
   * output is wriiten in the root file title.root
   */
void reaction::PlotPotentialEcm()
{
#ifdef root
  cout << "pot" << endl;
  TCanvas *pot = new TCanvas("pot");

  double Emin = -200.;
  double Emax = 200.;
  // axies
  TH2S *hist4 = new TH2S("hist4",title.c_str(),10,Emin,Emax
                        ,10,-15.,40.);
  hist4->GetXaxis()->SetTitle("E_{c.m.} [MeV]");
  hist4->GetYaxis()->SetTitle("Potential [MeV]");
  hist4->SetStats(kFALSE);
  hist4->Draw();
  int const Npoints =150;
  double ee[Npoints];
  double WWsurfaceLow[Npoints];
  double WWvolume[Npoints];
  double VVHF[Npoints];
  double VVsurfaceLow[Npoints];
  double VVvolume[Npoints];
  double VVHFsurface[Npoints];
  double Vtotal[Npoints];



  string file_name = title +".poten";
  ofstream pFile(file_name.c_str());
  for (int i=0;i<Npoints;i++)
    {
      double Elab = (Emax-Emin)/(double)Npoints*(double)i + Emin;
      double Ecm = energyLab2Cm(Elab);
      ee[i] = Ecm;

      complex< double > volumePot = Pot->volumeE( Ecm );
      complex< double > surfacePot = Pot->surfaceE( Ecm );

      WWvolume[i] = imag( volumePot );
      WWsurfaceLow[i] = imag( surfacePot );

      VVHF[i] = Pot->HartreeFock.Vvol;
      VVHFsurface[i] = Pot->HartreeFock.Vsur;

      VVvolume[i] = real( volumePot );
      VVsurfaceLow[i] = real( surfacePot );

      Vtotal[i] = VVHF[i] + VVvolume[i];

      pFile << Ecm-Efermi << " " << WWsurfaceLow[i]  
              <<   " " << WWvolume[i] << " " 
	    << VVHF[i] << " " << VVHFsurface[i] << " " << 
            VVsurfaceLow[i] <<
	" " << VVvolume[i] << endl;


      VVHF[i] /= 6.;
      Vtotal[i] /= 6.;




    }
  pFile.close();
  pFile.clear();

  TGraph *gWvolume = new TGraph(Npoints,ee,WWvolume);
  TGraph *gWsurfaceLow = new TGraph(Npoints,ee,WWsurfaceLow);
  TGraph *gVvolume = new TGraph(Npoints,ee,VVvolume);

  TGraph *gVsurfaceLow = new TGraph(Npoints,ee,VVsurfaceLow);
  TGraph *gVHF = new TGraph(Npoints,ee,VVHF);
  TGraph *gVHFsurface = new TGraph(Npoints,ee,VVHFsurface);
  TGraph *gVtotal = new TGraph(Npoints,ee,Vtotal);


  


  gWvolume->SetLineColor(2);
  gWsurfaceLow->SetLineColor(6);
  gWvolume->SetLineWidth(3);
  gWsurfaceLow->SetLineWidth(3);
  gVvolume->SetLineColor(2);
  gVsurfaceLow->SetLineColor(6);

  gVHF->SetLineColor(4);
  gVHFsurface->SetLineColor(5);
  gVtotal->SetLineColor(1);
  gVtotal->SetLineWidth(3);
  


  gWvolume->Draw("C");
  gWsurfaceLow->Draw("C");
  gVvolume->Draw("C");
  //gVisoscalerSurface->Draw("C");
  //gVisovectorSurface->Draw("C");
  gVsurfaceLow->Draw("C");


  gVHF->Draw("C");
  gVHFsurface->Draw("C");
  gVtotal->Draw("C");


  pot->Write();



  //------------------------------------------------------------
  cout << "wfunct" << endl;
  //plot wave functions
  TCanvas *wfunct = new TCanvas("wfunct");
  wfunct->Divide(2,2);
  TH2S *hist6 = new TH2S("hist6",title.c_str(),10,0.,35.
                        ,10,-.7,.7);

  hist6->GetXaxis()->SetTitle("radius [fm]");
  hist6->GetYaxis()->SetTitle("reduced Wave Function");
  hist6->SetStats(kFALSE);
  TLine *line0 = new TLine(0.,0.,35.,0.);
  for (int i=1;i<=4;i++)
    {
    wfunct->cd(i);
    hist6->Draw();
    line0->Draw();
    }

  TGraph *graphW[15];
  double rrr[BoundRspace->nWave];
  for (int i=0;i<BoundRspace->nWave;i++) 
    rrr[i] = BoundRspace->rStart + (double)i*BoundRspace->dx;




  //***here

 ofstream fw ("p12.dat"); //rjc


  NlevelTh = 0;
  int ifine = 1;
  int nmax[8] = {2,2,2,1,1,0,0,0};
  char cwave[8]={'s','p','d','f','g','h','i','j'};
  double Elower= -90.;
  double Eupper = (double)BoundRspace->LogDerMax;
  if (Zp == 1 && Eupper > 10.) Eupper = 10.;
  for (int l=0;l<8;l++)
    {
      for (int ii=-1;ii<2;ii+=2)
	{
          double j = (double)l + (double)ii*0.5;
          if (j < 0.) continue;
          Elower = -70.;           

          for (int n=0;n<=nmax[l];n++)
	    {
             if (BoundRspace->searchNonLoc(Elower,Eupper,j,l,ifine)==1.) 
                {
                  //Add on the exterior wave function
                  BoundRspace->exteriorWaveFunct(l);
		  //normalise the wavefunction
		  BoundRspace->normalizeWF();
		 ostringstream outstring;
                 outstring << n << cwave[l] <<"#frac{"<<(int)(2.*j)<<"}{2}";
		 string name = outstring.str();
                 cout << name << endl;

                 LevelTh[NlevelTh].j = j;
                 LevelTh[NlevelTh].l = l;
                 LevelTh[NlevelTh].N = n;
                 LevelTh[NlevelTh].color = l + 1;
                 LevelTh[NlevelTh].name = name;
                 LevelTh[NlevelTh].energy = Elower;
                 //LevelTh[NlevelTh].Rrms = scatter->Rrms;
                 graphW[NlevelTh] = new TGraph(BoundRspace->nWave,rrr,
                 BoundRspace->WaveArray);
                 graphW[NlevelTh]->SetLineColor(LevelTh[NlevelTh].color);
                 wfunct->cd(LevelTh[NlevelTh].l+1);
                 graphW[NlevelTh]->SetLineStyle(n+1);
                 graphW[NlevelTh]->Draw("C");
		 Elower = BoundRspace->Elower0 + 0.5;
                 NlevelTh++;


		 if (n ==1 && l == 1 && j == 0.5) //rjc
		   {
                     for (int i=0;i<BoundRspace->nWave;i++)
                        {
			  fw << rrr[i] << " " << BoundRspace->y[i] << endl;
                        }
		   }

		}
	    
	     else break;
	    }
	}
    }
 
  fw.close();  //rjc
  fw.clear();  //rjc

  int jj = 0;
  // sort levels in energy order
  for (int j=0; j<NlevelTh;j++)
    {
      double Emin = 10000000.;
      int imin=1000000;
      for (int i=0;i<NlevelTh;i++)
        {
	  if (LevelTh[i].energy < Emin)
	  {
	    Emin = LevelTh[i].energy;
            imin = i;
	  }
        }
      LevelThSort[jj] = LevelTh[imin];
      LevelTh[imin].energy = 10000000000.;
      jj++;
    }

  
  wfunct->Write();

  



  

  //-------------------------------------------------------
  cout << "lev" << endl;
  //now for the levels
  double const TheoryLeft = 1.1;
  double const TheoryRight = 1.7;
  double const ExpLeft = .3;
  double const ExpRight = .9;
  TCanvas *lev = new TCanvas("lev");
  TH2S *hist7 = new TH2S("hist7",title.c_str(),10,0.,2.
                        ,10,-70,5);

  hist7->GetYaxis()->SetTitle("Level Energy [MeV]");
  hist7->SetStats(kFALSE);
  hist7->Draw();
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextSize(.02);

  //draw fermi energy
  TLine *Lfermi = new TLine(ExpLeft,Efermi,TheoryRight,Efermi);
  Lfermi->SetLineWidth(3);
  Lfermi->SetLineStyle(2);
  Lfermi->Draw();

  TLine *Lf[NlevelTh];
    for (int i=0;i<NlevelTh;i++)
      {
      Lf[i] = new TLine(TheoryLeft,LevelThSort[i].energy,
                          TheoryRight,LevelThSort[i].energy);
      Lf[i]->SetLineColor(LevelThSort[i].color);
      Lf[i]->SetLineWidth(3);
      Lf[i]->Draw();
      tex.SetTextColor(LevelThSort[i].color);
      tex.DrawLatex(TheoryRight,LevelThSort[i].energy,
                          LevelThSort[i].name.c_str());
      }


  //plot experimental levels
  TLine line;
  line.SetLineWidth(3);
  if (Nlevel > 0) 
     {

     TLine *expLevel[Nlevel];
     for (int i=0;i<Nlevel;i++)
       {
         line.SetLineColor(Level[i].color);
         expLevel[i] = new TLine(ExpLeft,Level[i].Energy,
                                 ExpRight,Level[i].Energy);
         expLevel[i]->SetLineWidth(3);
         expLevel[i]->SetLineColor(Level[i].color);
         expLevel[i]->Draw();

         //connect to calculated value
         for (int j=0;j<NlevelTh;j++)
	   {
	     if (Level[i].l == LevelThSort[j].l &&
                 Level[i].j == LevelThSort[j].j &&
                 Level[i].N == LevelThSort[j].N)
	         {
		   line.DrawLine(ExpRight,Level[i].Energy,
		   			 TheoryLeft,LevelThSort[j].energy);
	         }
	   } 

       }
     }

  //write out in file level stuff
  string filename = title + ".level";
  ofstream LevelStuff(filename.c_str());
  LevelStuff << "N J  L       E     Rrms    occup   SpectF    Width color" << endl; 
    for (int i=0;i<NlevelTh;i++)
      {
	LevelStuff << LevelThSort[i].N << " " <<  LevelThSort[i].j << 
                   " " << LevelThSort[i].l 
		   << " "  << LevelThSort[i].energy << " " <<
	  LevelThSort[i].Rrms << " " << LevelThSort[i].Occupation << " " <<
	  LevelThSort[i].SpectFactor << " " << LevelThSort[i].Width << 
	  " "  << LevelThSort[i].ANC << " "   << LevelThSort[i].color << endl;
      }
  LevelStuff.close();
  LevelStuff.clear();


  lev->Write();

  //-----------------------
  cout << "react" << endl;
  TCanvas *Creact = new TCanvas("react");
  double ymax;
  double xmax;

  if (Zp > 0.) 
    {

      ymax = 1.;
      for (int i=0;i<NXdata;i++)
	{
	  if (Xdata[i].xsec > ymax) ymax = Xdata[i].xsec;
	}
      ymax *= 1.2;
     xmax = 200.;
    }
  else 
    {
      if (Z > 39) ymax = 10000.;
      else ymax = 4500.;
     xmax = 200.;
    }
  TH2S histXsec("histXsec",title.c_str(),10,0,xmax,10,0,ymax);

  histXsec.SetStats(kFALSE);
  histXsec.GetXaxis()->SetTitle("E_{lab} [MeV]");
  if (Zp > 0.) histXsec.GetYaxis()->SetTitle("#sigma_{react} [mb]");
  else histXsec.GetYaxis()->SetTitle("#sigma_{tot} [mb]");
  histXsec.Draw();


  double xxx[400];
  double yyy[400];
  double dxxx[400];
  double dyyy[400];



  filename = title+"r.dat";
  ofstream fdata(filename.c_str());

  //experimental
  fdata << NXdata << endl;
  if (NXdata > 0)
    {

     for (int i=0;i<NXdata;i++)
       {
	 fdata << Xdata[i].energyLab << " " <<
	   Xdata[i].xsec << " " << Xdata[i].sigma <<endl;
         xxx[i] = Xdata[i].energyLab;
         dxxx[i] = 0.;
         yyy[i] = Xdata[i].xsec;
         dyyy[i] = Xdata[i].sigma;
 
       }


     TGraphErrors *xsec_exp = new TGraphErrors(NXdata,xxx,yyy,dxxx,dyyy);
     xsec_exp->SetMarkerStyle(21);
     xsec_exp->SetMarkerColor(2);
     xsec_exp->Draw("P");
    }

  fdata << 143 << endl;

  //calculated
  for(int i=0;i<100;i++)
    {
      double Elab =((double)i+.5)/5. + .25;
      xxx[i] = Elab;
     double Ecm = energyLab2Cm(Elab);



       // integrate wavefunction and find phaseshifts
     ScatterRspace->getSmatrix(Ecm,Elab);

       yyy[i] = ScatterRspace->AbsorptionXsection();

       dxxx[i] = yyy[i] - ScatterRspace->statistical(ScatterRspace->konst,Ecm);

       fdata << Elab << " " << dxxx[i] <<endl;
    }


  TGraph xsec_th(100,xxx,yyy);
  xsec_th.Draw("C");
  xsec_th.SetLineColor(2);
  xsec_th.DrawGraph(100,xxx,dxxx);

  for (int i=0;i<43;i++)
    {
      double Elab = 20. + (double)i*180./32.;
      double Ecm = energyLab2Cm(Elab);

      xxx[i] = Elab;


      ScatterRspace->getSmatrix(Ecm,Elab);
      yyy[i] = ScatterRspace->AbsorptionXsection();

      fdata << Elab << " " << yyy[i] << endl;
    }


  TGraph *gxsec2 = new TGraph(33,xxx,yyy);
  gxsec2->SetLineColor(2);
  gxsec2->Draw("C");

  cout << "totxsec" << endl;
  //total cross sections for neutrons
  if (Zp == 0. && NTotXdata > 0)
    {
      fdata << NTotXdata << endl;
     for (int i=0;i<NTotXdata;i++)
       {
	 fdata << TotXdata[i].energyLab << " " <<
	   TotXdata[i].xsec << " " << TotXdata[i].sigma <<endl;
         xxx[i] = TotXdata[i].energyLab;
         dxxx[i] = 0.;
         yyy[i] = TotXdata[i].xsec;
         dyyy[i] = TotXdata[i].sigma;
 
       }


     TGraphErrors *txsec_exp = new TGraphErrors(NTotXdata,xxx,yyy,dxxx,dyyy);
     txsec_exp->SetMarkerStyle(22);
     txsec_exp->SetMarkerColor(4);
     txsec_exp->Draw("P");
    }

  cout << " ll " << endl;

  if (Zp == 0.)
    {
     //calculated
      fdata << 100 << endl;
     for(int i=0;i<100;i++)
       {
         double Elab =(double)(i+1)*2.;

         double Ecm = energyLab2Cm(Elab);
         xxx[i] = Elab;



         // integrate wavefunction and find phaseshifts
         ScatterRspace->getSmatrix(Ecm,Elab);

         yyy[i] = ScatterRspace->TotXsection();

         fdata << Elab << " " << yyy[i] << endl;
         //cout << Elab << " " << yyy[i] << endl;
       }


    TGraph *txsec_th = new TGraph(100,xxx,yyy);
    txsec_th->SetLineColor(4);
    txsec_th->Draw("C");
    }

  fdata.close();
  fdata.clear();



  Creact->Write();

 
    
#endif

}
//**************************************************************************
  /**
   * opens root file to save spectra in
   */
void reaction::OpenRootFile()
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
   * closes the root file
   */
void reaction::CloseRootFile()
{
#ifdef root
  f->Write();
  f->Close();
#endif
}

