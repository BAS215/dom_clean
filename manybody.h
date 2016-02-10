//
#ifndef _MANYBODY_H_
#define _MANYBODY_H_

#include <gsl/gsl_sf_bessel.h>
#include <list>
#include <utility>

#include "eigen.h"
#include "read_parameters.h"
#include "hamiltonian.h"

typedef std::pair< double, std::vector<double > > eigen_t;


cmatrix_t propagator( const std::vector<double> &rmesh, double E, int l, 
                      double j, const NuclearParameters &Nu, pot &U );

std::vector< double >
bound_eigvals( const matrix_t &Ham );

std::vector< double >
real_eigvals( const matrix_t &Ham ); 

bool 
compare_N( std::pair< double, std::vector<double> > pair1, 
           std::pair< double, std::vector<double> > pair2 );

std::vector< eigen_t >
bound_eigvecs( const matrix_t &Ham );

std::vector< eigen_t >
real_eigvecs( const matrix_t &Ham ); 

std::vector<eigen_t>
slfcn_eig( const std::vector<double> &rmesh, double Estart, int nb, int l, 
           double j, const NuclearParameters &Nu, pot &U );

double
find_level( const std::vector<double> &rmesh, double Estart, int N, int l, 
            double j, double tol, const NuclearParameters &Nu, pot &U );

eigen_t
find_boundstate( const std::vector<double> &rmesh, double QPE, int N, int l, 
                 double j, double tol, const NuclearParameters &Nu, pot &U );

std::vector< double > 
normalize( const std::vector< double > &rmesh, 
           const std::vector<double> &eigenvector );

double 
sfactor( const std::vector<double> &rmesh, double QPE, int l, double xj, 
         const std::vector<double> &QPF, pot &U );

double rms_radius( const std::vector<double> &rmesh, 
                   const std::vector<double> &QPF );

double occupation( const std::vector<double> &rmesh, const matrix_t &d_mtx,
                   const std::vector<double> &QPF );

double 
width( const std::vector<double> &rmesh, double QPE, int L, double J,
       const std::vector< double > &QPF, pot &U ); 

#endif // _MANYBODY_H_
