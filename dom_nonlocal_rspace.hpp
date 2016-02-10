/*****************************************************

        dom_nonlocal_r class

     Declarations for non-local dom in r-space
     and transformation to k-space


 *****************************************************/
#ifndef dom_nonlocal_rspace_hpp
#define dom_nonlocal_rspace_hpp
#include "pot.h"
#include "read_parameters.h"
#include "legendre.h"
#include "boundRspace.h"
#include "meshes.h"
#include "numerical.h"
#include "io.h"
#include "FourierBessel.hpp"
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <Eigen/LU>
#include <Eigen/Dense>
#include "Gauss_q.h"
#include "types.h"


using namespace std;
using namespace Eigen;
//
// Some constants of the DOM that do not seem to change
// Decided to put them here since they should not change in
// scattering process.

  //Potential specifications
   int const type = 1; // 1 is P.B. form (average), 0 is V.N. form
   int const mvolume = 4;  //power of the term in volume potential
   int const AsyVolume = 1;

//Neutron
    double const Zp =0.;// 1.;
   double const tz = -.5;//.5; 
//Proton
 //  double const Zp = 1.;
 //  double const tz = .5; 

   complex<double> const zero = complex<double> (0.,0.);
 class dom_nonlocal_r
 {
  private:
      double Energy_cm;
      int n;
      int l;
      double J;
      vector<double> rmesh;
      vector<double> drmesh;
  public:
      dom_nonlocal_r( double Ecm, 
                      int radialQuantumNumber,
                      int OrbitalAngularL, 
                      double TotalAngularJ, 
                      vector<double> &r,
                      vector<double> &dr);
 

     vector<double> dom_r_space_wavefunction(vector<double> rmesh );
  //   MatrixXcd dom_r_space(vector<double> &wf);
     MatrixXcd dom_r_space(cvector_t &Vdomlocal, vector<double> &wf);

     MatrixXcd dom_k_space( vector<double> &k, vector<double> &wf_k);

   //  MatrixXcd dom_r_space();
     MatrixXcd dom_r_space( cvector_t &Vdomlocal );

     MatrixXcd dom_k_space( vector<double> &k);
     int M_initial;
};

#endif
