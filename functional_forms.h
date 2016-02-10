#ifndef _FUNCTIONAL_FORMS_H_
#define _FUNCTIONAL_FORMS_H_

#include "types.h"
#include "Legendre.h"
#include "constants.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <cmath>
#include <utility>

double w( int, double, double, double, double, double );

double lorentzian( double x, double x0, double width, double norm ); 

double gaussian( double x, double x0, double sigma, double norm ); 

double woods_saxon( double, double, double );
double der_woods_saxon( double, double, double );
double der_der_woods_saxon( double, double, double );

double q_woods_saxon( double q, double R0, double a0 ); 

double 
volume_nl_form_PB1( double r1, double r2, double x, double R0, double a ); 

double 
volume_nl_form_PB2( double r1, double r2, double R0, double a ); 

double
surface_nl_form_PB2( double r1, double r2, double R0, double a ); 

double 
volume_nl_form_VN( double r1, double r2, double R0, double a ); 

double
surface_nl_form_VN( double r1, double r2, double R0, double a ); 

double gaussian_nl( double r1, double r2, double x, double beta );

matrix_t 
volume_nl_PB1( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta, const LegendrePoly &Pl ); 
matrix_t 
q_volume_nl_PB1( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta, const LegendrePoly &Pl ); 

matrix_t 
volume_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
                double a0, double beta ); 
matrix_t 
q_volume_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
                 double a0, double beta ); 
matrix_t 
surface_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
                double a0, double beta ); 

matrix_t 
volume_nl_VN( const std::vector<double> &rmesh, int l, double R0, 
              double a0, double beta ); 
matrix_t
surface_nl_VN( const std::vector<double> &rmesh, int l, double R0,
               double a0, double beta ); 
matrix_t 
q_form_bessel( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta ); 

double step( double diff ); 
double energyLab2Cm(double Elab, double A );
double energyCm2Lab(double Ecm, double A );

double WaveNumber( double energyCM, double A, double EminRel ); 
double Rel_Factor( double energyCM, double EminRel );

// Dispersion Integrals
double delta_w2_hole( double, double, double, double, double );
double delta_w2_particle( double, double, double, double, double );
double delta_w4_hole( double, double, double, double, double );
double delta_w4_particle( double, double, double, double, double );
double delta_asymmetric_hole( double, double, double, double );
double delta_asymmetric_particle( double, double, double, double );
double delta_w2( double, double, double, double, double );
double delta_w4( double, double, double, double, double );

// Derivatives of Dispersion Integrals
double der_delta_w2_hole( double, double, double, double, double );
double der_delta_w2_particle( double, double, double, double, double );
double der_delta_w4_hole( double, double, double, double, double );
double der_delta_w4_particle( double, double, double, double, double );
double der_delta_asymmetric_hole( double, double, double, double );
double der_delta_asymmetric_particle( double, double, double, double );


#endif
