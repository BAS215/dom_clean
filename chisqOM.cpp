// to run this program 
// chisq title file jdata X  
// title = name of inp file (without extension) containing OM parameter
// file =  name of reaction
// jdata = data set in file.data to fit
//X = x to fit or x=a to make calculation with given input parameter




#include "fitOM.h"
#include "compoundException.h"

int main(int Narg, char* pszArg[])
{

  if (Narg > 5) return 0;
  string title(pszArg[1]);
  string reactionTitle(pszArg[2]);
  string fileName = title+".inp";
  int jdata = atoi(pszArg[3]);
  string toFit(pszArg[4]);  

  string fileO = reactionTitle+".vvv";
  ofstream fileOut(fileO.c_str());
  fileOut.precision(3);
  fileOut << fixed;

  fileO = reactionTitle+"A.grid";
  ofstream fileGrid(fileO.c_str());
  fileGrid.precision(3);

  fileO = reactionTitle+"A.x";
  ofstream fileXsec(fileO.c_str());
  fileXsec.precision(3);
  fileXsec<<fixed;

  fileO = reactionTitle+"A.min";
  ofstream fileMin(fileO.c_str());
  fileMin.precision(3);


  //read in input file to cound how many fitted parameters
  ifstream file(fileName.c_str());
  if (file.fail())
    {
      cout << "could not open file " << fileName << endl;
      abort();
    }
  //char line[80];
  //skip two lines
  //file.getline(line,80);
  string variable;
  double parameter,varied,squared,scaling,pMin,pMax;
  int Ndim = 0;
  for (int i=0;i<fitOM::TotPara;i++)
    {
      file >> parameter >> varied >> squared>>scaling >> pMin>> 
              pMax >>variable;
      if (varied) Ndim++;
    }  
  

  file.close();
  file.clear();
  fitOM *Fit;

  //  for (int id =0;id<=jdata;id++)

  //    {

       int id = jdata;
      cout << "fitting data set " << id << endl;
      Fit = new fitOM (&title,reactionTitle,id,Ndim,0,1);

     //end of reaction data
     if (Fit->React[0]->data[0].nX < 0) 
       {
	 delete Fit;
	 return 0;
       }

     // if there are no xsec data , then skip
     if (Fit->React[0]->data[0].nX == 0) 
       {
         cout << "no xsec data" << endl;
         delete Fit;
         return 0;
       }

     //check if number of fit parametersed changed
     if (Ndim != Fit->mm)
       {
	 Ndim = Fit->mm;
       }
     delete Fit;
     Fit = new fitOM(&title,reactionTitle,id,Ndim,0,1);
   



     //grid search
     double chiMin=0.;
     if (toFit == "G" || toFit ==  "g")
       {
         chiMin = Fit->gridFit();
         cout << " min chisq from grid = " << chiMin << endl;
         fileGrid << id << " " << Fit->React[0]->data[0].energyLab << " " 
                  <<chiMin << " " ;
         for (int i=0;i<Fit->TotPara;i++) fileGrid << Fit->paraBest[i] << " ";
	 fileGrid <<  Fit->React[0]->data[0].name.substr(0,9) << endl;

       }

  
     double * para;
     para = new double [Ndim];

     for (int i=0;i<Ndim;i++)
       {
         int j= Fit->map2[i];
         para[i] = Fit->allPara[j]*Fit->scaling[j];
       }

      double* xi[Ndim];
     for (int i=0;i<Fit->ND;i++) 
       {
         xi[i] = new double [Ndim];
         for (int j=0;j<Fit->ND;j++)
	   {
	     if (i == j) xi[i][j] = 1.;
	     else xi[i][j] = 0.;
	   }
       }
     double const ftol = 1.e-2;



     //string toFit("X");
     cout << toFit << endl;
       if (toFit == "X" || toFit ==  "x" || toFit == "G" || toFit == "g")
       {
         chiMin=Fit->powell(para,xi,ftol);
         cout << Fit->React[0]->chidif << endl;
       }

     cout << "minimum chisq= " << chiMin << endl;

     Fit->SavePara(para);

     Fit->WriteFit(chiMin);

     for (int i=0;i<Fit->Nreact;i++)
       {
	 Fit->React[i]->OpenRootFile();
	 cout << "PlotFit" << endl;
         try
	   {
	    Fit->React[i]->PlotFit();
            cout<< "absorption xsec = " << Fit->React[i]->scatter->AbsorptionXsection() << endl;
	    cout << "kwave = " << Fit->React[i]->scatter->Kwave << endl;
	   }
         catch (compoundException)
	   {
	     cout << "compound exception" << endl;
	   }
         //Fit->React[i]->plotSmatrix();
	 Fit->React[i]->CloseRootFile();

       }

     cout << " for Ecm = "<< Fit->React[0]->data[0].energyCM << " MeV" << endl;
     cout << " for Elab = "<< Fit->React[0]->data[0].energyLab << " MeV" << endl;
     //     if (Fit->React[0]->Zp == 0) 
     //  cout << " predicted= " << Fit->React[0]->scatter->TotXsection() << 
     //  "exp = " << Fit->React[0]->TotXdata[0].xsec << " mb "<< endl;

     Fit->React[0]->scatter->VIntegrals();
     cout << " JReal= " << Fit->React[0]->scatter->JReal 
          << " JImag= " << Fit->React[0]->scatter->JImag 
          << " RrmsReal= " << Fit->React[0]->scatter->RrmsReal 
          << " RrmsImag= " << Fit->React[0]->scatter->RrmsImag << endl;

     fileOut<< setw(5) << Fit->React[0]->data[0].energyLab << " " 
            << setw(6) << -Fit->React[0]->scatter->JReal  << " " 
	    << setw(6) <<-Fit->React[0]->scatter->JImag << " " 
            << setw(5) << Fit->React[0]->scatter->RrmsReal  << " " 
	    << setw(5) << Fit->React[0]->scatter->RrmsImag  << " " ;
     if (Fit->React[0]->data[0].nA == 0)fileOut << setw(6) << 0. << " " 
                                                << setw(5) << 0. << " ";
     else fileOut << setw(6) << -Fit->React[0]->scatter->JSO << " " 
		  << setw(5) << Fit->React[0]->scatter->RrmsSO << " " ;
     fileOut << setw(5) << chiMin << " " << 
       Fit->React[0]->data[0].name.substr(0,9) <<endl;

     /*
     Fit->React[0].loadOM();
     Fit->React[0]->scatter->statistical(scatter->konst,
         Fit->React[0]->data[0].energyEcm);
     */
     fileXsec << setw(3) << id << " " 
             << setw(6) << Fit->React[0]->data[0].energyLab << " "
	     << setw(6) << Fit->React[0]->scatter->AbsorptionXsection()<<  " "
	     << setw(6) << Fit->React[0]->scatter->sigmaAbsorption<< " " 
	     << setw(6) << Fit->React[0]->scatter->sigmaCompoundElastic<< " " ;
           
     if (Fit->React[0]->Zp == 0) fileXsec << setw(6) << 
       Fit->React[0]->scatter->TotXsection() << " " ;

     fileXsec   <<  Fit->React[0]->data[0].name.substr(0,9) << endl;


     fileMin << id << " " << Fit->React[0]->data[0].energyLab << " "
	     << chiMin << " "; 
     for (int i=0;i<Fit->TotPara;i++) fileMin << Fit->allPara[i] << " " ;
     fileMin<<  Fit->React[0]->data[0].name.substr(0,9) << endl;
     delete Fit;




     // }
     fileOut.close();
     fileOut.clear();

     fileGrid.close();
     fileGrid.clear();

     fileXsec.close();
     fileXsec.clear();

     fileMin.close();
     fileMin.clear();


  return 0;
}

