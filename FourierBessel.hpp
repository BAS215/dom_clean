/*
    Simple class that contains the 
   spherical bessel function values at 
   all mesh points k-r
    Members of the class:
     l := Orbital Angular Momentum Quantum number
     q := Either momentum or position radial component

    j_l_kr()
      Input:
         p:= Canonical conjugated variable of q,
             e.g. q is momentum, p is position
      Output:
        j_l(qp) := spherical bessel function at 
                   q times p
 */
 #ifndef FourierBessel_hpp
 #define FourierBessel_hpp
 #include <iostream>
 #include <cmath>
 #include <vector>
 #include <boost/math/special_functions.hpp>
 #include <Eigen/LU>
 #include <Eigen/Dense>
  using namespace std;
  using namespace boost::math;
  using Eigen::MatrixXd;
  using Eigen::MatrixXcd;

 class FourierBesselContainer
 {
   private:
       int l;
       vector<double> q;
   public:
     FourierBesselContainer(int OrbitalAngularL,vector<double> &QMesh)
     { l= OrbitalAngularL; q = QMesh; }
     vector<double> j_l_kr(double p);
     MatrixXd all_j_l_kr(vector<double> &p);
 };
#endif
