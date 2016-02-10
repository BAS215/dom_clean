/***************************************************

   ReducibleSelfEnergy Class declarations

 ***************************************************/
 #ifndef RedSelfEnergy_hpp
 #define RedSelfEnergy_hpp
 #include "constants.hpp"
 #include "dom_nonlocal_rspace.hpp"
 #include <iostream>
 #include <boost/math/special_functions.hpp>
 #include <Eigen/LU>
 #include <Eigen/Dense>
 using namespace std;
 using Eigen::MatrixXd;
 using Eigen::MatrixXcd;
 double const hbc2=hbarc*hbarc;
 class ReducibleSelfEnergy
 {
  private:
     int n;
     int l;
     double j;
     double Ecm;
     double Mu;
     vector<double> k;
     vector<double> dk;
     vector<double> Go;
  public:
    ReducibleSelfEnergy(int radialQuantumN,
                         int OrbAngularL,
                         double TotAngularJ,
                         double EnergyCM,
			 double ReducedMass,
                         vector<double> &wvVector,
                         vector<double> &dwvVector,
                         vector<double> &FreeProp);

   void FindRMatrix( MatrixXcd &R, MatrixXcd &cR, vector<double> &wf);
   void FindRMatrix( MatrixXcd &R, MatrixXcd &cR);

   void FindRedSigma(MatrixXcd &SelfEnergy, MatrixXcd &cSelfEnergy, 
                     vector<double> &wf);

   void FindRedSigma(MatrixXcd &SelfEnergy, MatrixXcd &cSelfEnergy); 


   void SelfEnergyAtPole(complex<double> &S_lj,complex<double> &cS_lj);

   //added togenerate bound state wf to calculate S(E) using Sp_rr in spectralFunction.cpp 
   vector<double> FindRedSigma( MatrixXcd &R, MatrixXcd &cR , vector<double> r,vector<double> dr);

 };
#endif
