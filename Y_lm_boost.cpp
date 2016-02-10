#include </usr/local/boost_1_48_0/boost/math/special_functions/spherical_harmonic.hpp>
#include <gsl/gsl_sf_legendre.h>
#include "legendre.h"
#include "boundRspace.h"
#include <gsl/gsl_complex.h>
#include "pot.h"
#include <cmath>
int main() {

 int l = 5;
 int m = 1;
 //int YYlm_m = (int) Ylm_m;
 double theta = M_PI/4.;
 double phi = M_PI/13.;
 double angle = m * phi; 
// std::complex<double> exponent(cos(angle),sin(angle));
// complex<double> Polym(gsl_sf_legendre_sphPlm(l,m,theta),0.0);
// double Polym(gsl_sf_legendre_sphPlm(l,m,theta));
// std::cout<< exponent * Polym <<std::endl;
// std::complex<double> BUSS( boost::math::spherical_harmonic(l,m, theta, phi));
// std::cout<< BUSS <<std::endl;
}

