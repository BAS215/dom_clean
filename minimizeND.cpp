#include "minimizeND.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

//************************************************************
minimizeND::minimizeND(int ND0)
{
  ND = ND0;
  pcom = new double [ND];
  xicom = new double [ND];
} 
//*************************************************************
minimizeND::~minimizeND()
{
  delete [] pcom;
  delete [] xicom;
}
//*************************************************************
//used by linmin. Provides a 1D function for use by the base class
//minimize1D
double minimizeND::funct1D(double x)
{
  double xt[ND];
  for (int i=0;i<ND;i++) xt[i] = pcom[i] + x*xicom[i];
  return functND(xt);
}
//**************************************************************
// minimize a function along a direction specified by the vector xi[],
// starting at the initial point p[]. The minimized function value is 
// returned. Used by function powell
double minimizeND::linmin(double*p,double*xi)
{
  for (int i=0;i<ND;i++)
    {
      pcom[i] = p[i];
      xicom[i] = xi[i];
    }
  //initial guess for brackets
  double xb[3],fb[3];
  xb[0] = 0.;
  xb[1] = 1.;
  bracket(xb,fb);
  //cout << "brackets = " << xb[0] << " " << xb[1] << " " << xb[2] << endl;
  //cout << "f= " << fb[0] << " " << fb[1] << " " << fb[2] << endl;
  double const tolerance = 1.e-2;
  double fret = Brent(xb,tolerance);
  for (int i=0;i<ND;i++)
    {
      xi[i] *= xb[0];
      p[i] += xi[i];
    }
  return fret;
}
//***************************
//multidimension minimization using the powell method.
//this was translated from "numerical recipes in Fortran 77" 
// by Press, Teukolsky, Vetterling, Flannery , 2nd edition page 411
double minimizeND::powell(double *p,double **xi,double const &ftol)
{
  const double Tiny=1.e-25;
  const int iterMax = 30;
  double fret = functND(p);
  double pt[ND],ptt[ND];
  double xit[ND];
  //cout << p[0] << " " << p[1] << endl;
  //save initial point
  for (int i=0;i<ND;i++) pt[i] = p[i]; //save initial point
  int iter = 0;

  std::ofstream chisqfiler( "Output_test/chisq_rec.out" );
  int mmm = 0;
  for (;;)
    {
      cout << "********iteration= " << iter << " funct= " << fret << endl;
     iter++;
     double fp = fret;
     int ibig = 0;
     double del= 0.;   //del will record the largest function decrease
     for(int i=0;i<ND;i++)
       {
         cout << "******para i = " << i << endl;
	 for (int j=0;j<ND;j++) xit[j] = xi[j][i];
         double fptt = fret;
         fret = linmin(p,xit);   // minimize along the direction xit
         //cout << fret << " " << p[0] << " " << p[1] << endl;
	 if (fptt-fret > del)  // check to see if this was the largest decrease
	   {
	     del = fptt - fret;
	     ibig = i;
	   }
	   mmm++;	
           chisqfiler <<mmm<<" "<< fret << std::endl;
       }

     if (2.*(fp-fret) <= ftol*(abs(fp)+abs(fret))+Tiny) return fret;
     if (iter == iterMax) 
       {
	 cout << "powell excedding maximum iterations" << endl;
         return fret;
       }
     //construct the extrapolated point and the average direction moved.
     //save the old starting place
     for (int j=0;j<ND;j++)
       {
	 ptt[j] = 2.*p[j]-pt[j];
	 xit[j] = p[j] - pt[j];
	 pt[j] = p[j];
       }
     double fptt = functND(ptt); // function value at extrapolated point
     if (fptt >= fp) continue;  // one reason not to use the new direction

     //The following lines were in the original code, 
     //but I am now sure what they are supposed to do
     //they certainly stops us going to a better minimumsome times

     //double t = 2.*(fp-2.*fret+fptt)*pow(fp-fret-del,2)-del*pow(fp-fptt,2);
     //if (t >= 0.) continue; // another reason not to use the new direction


     fret = linmin(p,xit); // move to the minimum of the new direction
     //save new direction
     for (int j=0;j<ND;j++)
       {
	 xi[j][ibig] = xi[j][ND-1];
         xi[j][ND-1] = xit[j];
       }
    }
   chisqfiler.close();
}
