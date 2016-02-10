/******************************************************************
   Spectral Function programme
            --------------

  File: spectralFunction.cpp
  Programmer:Helber Dussan
  Contact:hdussan AT physics.wustl.edu

       Particle Spectral Function Class
 
  Has option to use either reduced mass
  or nucleon mass.
 ______________________________________________________
  void get_ParticleSpectralFunction(double Mu,int n, int l, double j)
 
   High level function to calculate the 
  Particle Spectral Density
     INPUT PARAMETRES
   Mu : Reduced Mass in MeV
   n  : Princ.Quantum Number
   l   : Orbital angular momentum
   j   : Total angular momentum

  - It does not keep track of each contribution
  - Using spectralFunction()
 ______________________________________________________
  void get_ParticleSpectralFunction2(double Mu,int n, int l, double j)
  
  High level function to calculate the 
  Particle Spectral Density
     INPUT PARAMETRES
   Mu : Reduced Mass in MeV
   n  : Princ.Quantum Number
   l   : Orbital angular momentum
   j   : Total angular momentum

  - It keeps track of each contribution
  - Using spectralFunction2()

 ______________________________________________________
 complex<double> spectralFunction(double Energy_cm,
                                  int n, int l, double j)
      INPUT PARAMETRES
   Ecm : Centre of Mass Energy
   n   : Principal quantum number
   l   : Orbital angular momentum
   j   : Total  angular momentum
      OUTPUT
  Particle Spectral Function

 ______________________________________________________

 vector< complex<double> > spectralFunction2(double Energy_cm, 
                                             int n, int l, double j)

      INPUT PARAMETRES
   Ecm : Centre of Mass Energy
   n   : Principal quantum number
   l   : Orbital angular momentum
   j   : Total  angular momentum
      OUTPUT
  Particle Spectral Function
  Each contribution is contained in the components of the vector
   Total Spectral function is contained in the 5th index.
 ---------------------------------------------------------

 MatrixXcd Sp_rr() : Full Spectral function in position space.
                     Functions Sp0_rr, Sp1_rr, Sp2_rr, etc. 
                     calculate separately each contribution, see
                     equation 0.30 in notes spectralfunction.pdf 
    Units:
      MeV/fm *fm2/fm2 =  MeV fm2 / fm3
                      =  MeV / (hbc2) MeV2 fm3 
                      =    1 / (hbc2 MeV fm3 )
**************************************************************/
#ifndef spectralFunction_hpp
#define spectralFunction_hpp
#include "constants.hpp"
#include "kMesh.hpp"
#include "RedSelfEnergy.hpp"
#include "FourierBessel.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

//complex<double> const zero = complex<double>(0.0,0.0);
class p_spectralFunction
{
  private:
    double Mu;
    int n;
    int l;
    double j;
  public:

    p_spectralFunction(double redMass,
                       int PrincQuantumN,
                       int angularL, 
                       double angularJ)
    {
      Mu =  redMass; 
      n  =  PrincQuantumN;
      l  =  angularL;   
      j =   angularJ;
    }


    void get_ParticleSpectralFunction();
   
   
    void get_ParticleSpectralFunction2();

   void savespectral(double Energy_cm);

   vector< complex<double> > depletion(vector<double> &k,
				      vector<double> &dk,
				      vector<double> &Go,
				      MatrixXcd &dsigma_plus,
				      MatrixXcd &dsigma_menos,
                                      vector<double> &wf);
   
    complex<double> spectralFunction(double Energy_cm);
   
    vector< complex<double> > spectralFunction2(double Energy_cm);


   //Transf. to r-space. Needs to be tested //
    MatrixXcd Sp0_rr(double ko,MatrixXd &J_kr);

 
    MatrixXcd Sp1_rr(int rsize,
                     vector<double> &dk,
                     vector<double> &k,
                     vector<double> &G0,
		     MatrixXcd &sigmakkMenos,
                     MatrixXd &J_kr);

    MatrixXcd Sp2_rr(int rsize,
                     vector<double> &dk,
                     vector<double> &k,
                     vector<double> &G0,
		     MatrixXcd &sigmakkPlus,
                     MatrixXd &J_kr);

    vector< complex<double> > Sp0wave_r(int rsize,
                     vector<double> &dk,
                     vector<double> &k,
                     vector<double> &G0,
		     MatrixXcd &sigmakkPlus,
                     MatrixXd &J_kr);

    vector< complex<double> > Sp1wave_r(int rsize,
                     vector<double> &dk,
                     vector<double> &k,
                     vector<double> &G0,
		     MatrixXcd &sigmakkPlus,
                     MatrixXd &J_kr);

    MatrixXcd Sp3_rr(int rsize,
                     vector<double> &dk,
                     vector<double> &k,
                     vector<double> &G0,
		     MatrixXcd &sigmakkPlus,
                     MatrixXd &J_kr);

    MatrixXcd Sp4_rr(double ko,
		     MatrixXcd &sigmakkMenos,
		     MatrixXd &J_kr);

    MatrixXcd Sp_rr(double Energy_cm,vector<double> r, vector<double> dr);
  //  MatrixXcd Sp_rr(double Energy_cm,vector<double> r);

};
#endif
