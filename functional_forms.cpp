/* Functions for calculating the DOM potential
*/

#include "functional_forms.h"

// 'w' calculates the forms used in imaginary potential
double w( int n, double E, double Ef, double a, double b, double c ) {

    double aemc = std::abs( E - Ef ) - c;
    return a * step( aemc ) * std::pow(aemc, n ) 
             / ( std::pow(aemc,n) + std::pow(b,n));
}

// returns a lorentzian normalized to 'norm'
double lorentzian( double x, double x0, double width, double norm ) {

    return norm * width / M_PI / ( std::pow( x - x0, 2 ) + std::pow( width, 2 ) );

}

// returns a Gaussian normalized to 'norm'
double gaussian( double x, double x0, double sigma, double norm ) {

    return norm * std::exp( - std::pow( ( x - x0 ) / sigma, 2 ) / 2 )
                / ( sigma * std::sqrt( 2 * M_PI ) );
}

// Form Factors. R0 is the parameter defining the radial 
// extent of the nucleus

// Woods-Saxon
double woods_saxon( double r, double R0, double a0 ) {
    
    return 1 / ( 1. + std::exp(( r - R0 ) / a0 ));
}

// Derivative of a Woods-Saxon
double der_woods_saxon( double r, double R0, double a0 ) {

    return 1 / ( - 2.* a0 * ( 1 + std::cosh( ( r - R0 ) / a0 )));
}

// 2nd Derivative of a Woods-Saxon
double der_der_woods_saxon( double r, double R0, double a0 ) {

    double theta = ( r - R0 ) / a0;
    return std::tanh( theta / 2 ) / ( 1 + std::cosh( theta ) ) / std::pow(a0,2);
}

double q_woods_saxon( double q, double R0, double a0 ) {

    double term1;
    if ( q == 0 ) {
        term1 = R0 * ( R0 * R0 + std::pow( a0 * M_PI, 2 ) ) 
              / ( 6 * M_PI * M_PI );
    }
    else {
        term1 = a0 / ( 2 * M_PI * q ) 
              * ( M_PI * a0 * std::sin( q * R0 ) * std::cosh( M_PI * a0 * q ) 
              / std::pow( std::sinh( M_PI * a0 * q ), 2 )
              - R0 * std::cos( q * R0 ) / std::sinh( M_PI * a0 * q ) );
    }

    double term2 = 0;
    int nmax = 3;
    for ( int n = 1; n < nmax; ++n ) {
    
        term2 -= std::pow( a0, 3 ) / (M_PI * M_PI )
               * std::pow( -1.0, n ) * n * std::exp( -n * R0 / a0 )
               / std::pow( n * n + std::pow( q * a0, 2 ), 2 );
    }

    return term1 + term2;
}

// --- Functions for defining nonlocal potential --- //

// Form used by Perey and Buck in Nuclear Physics 32 (1962) 353-380
// for the volume part of the potential. 
// The function below is just U(0.5*|\vec{r1} + \vec{r2}|). It 
// does not include the Gaussian factor
// x is cos( \theta )
double 
volume_nl_form_PB1( double r1, double r2, double x, double R0, double a ) {

    double rs = std::sqrt( std::pow( r1, 2 ) + std::pow( r2, 2 ) 
                           + 2 * r1 * r2 * x );

    return woods_saxon( rs / 2.0, R0, a );
}

// Form used by Perey and Buck in Nuclear Physics 32 (1962) 353-380
// for the volume part of the potential. 
// The function below is similar to the function above, but using
// the approximation 0.5*|\vec{r1} + \vec{r2}| = 0.5 * ( r1 + r2 ),
// so there is no angle dependence. 
double 
volume_nl_form_PB2( double r1, double r2, double R0, double a ) {

    double rs = 0.5 * ( r1 + r2 );

    return woods_saxon( rs, R0, a );
}

// Form used by Perey and Buck in Nuclear Physics 32 (1962) 353-380
// for the surface part of the potential. 
// The approximation 0.5*|\vec{r1} + \vec{r2}| = 0.5 * ( r1 + r2 ) is 
// again used, so there is no angle dependence. 
double
surface_nl_form_PB2( double r1, double r2, double R0, double a ) {

    double rs = 0.5 * ( r1 + r2 );

    return der_woods_saxon( rs, R0, a );
}

// Form proposed by Dmitri Van Neck. Consists of the square root
// of two Woods-Saxon form factors multiplied together.
// There is no angle dependence
double 
volume_nl_form_VN( double r1, double r2, double R0, double a ) {

    return std::sqrt( woods_saxon( r1, R0, a ) * woods_saxon( r2, R0, a ) );
}

// Possible form for the surface form factor, using a similar
// form as the previous function (volume_nl_form_VN)
double
surface_nl_form_VN( double r1, double r2, double R0, double a ) {

    return std::sqrt( der_woods_saxon( r1, R0, a ) * 
                      der_woods_saxon( r2, R0, a ) );
}

// This is the Gaussian factor, which determines the degree of nonlocality
// x is cos( theta )
double gaussian_nl( double r1, double r2, double x, double beta ) {

    double rv = std::sqrt( std::pow( r1, 2 ) + std::pow( r2, 2 ) 
                           - 2 * r1 * r2 * x );
    
    return std::exp( - std::pow( rv, 2 ) / std::pow( beta, 2 ) ) 
           / ( std::pow( M_PI, 1.5 ) * std::pow( beta, 3 ) );
}

// Perey-Buck representation of nonlocality
// Uses the form WS( \vec{r) + \vec{r}^\prime ) * Gaussian
// A numerical integration is required, since the Woods-Saxon
// term also contains angle dependence
matrix_t 
volume_nl_PB1( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta, const LegendrePoly &Pl ) {

    // Prepare Angular Integration
    int npart;
    std::vector<double> xil;
    std::vector<double> xiu;

    if( l == 0 ) {
        npart = 1;
        xil.push_back( -1 );
        xiu.push_back( 1 );
    }
    else if ( l == 1 ) {
        npart = 2;
        xil.push_back( -1 );
        xiu.push_back( 0 );
        xil.push_back( 0 );
        xiu.push_back( 1 );
    }
    else if ( l == 2 ) {
        npart = 3;
        xil.push_back( -1 );
        xiu.push_back( -1 / std::sqrt( 3 ) );
        xil.push_back( -1 / std::sqrt( 3 ) );
        xiu.push_back( 1 / std::sqrt( 3 ) );
        xil.push_back( 1 / std::sqrt( 3 ) );
        xiu.push_back( 1 );

    }
    else {
        npart = 2;
        xil.push_back( -1 );
        xiu.push_back( 0 );
        xil.push_back( 0 );
        xiu.push_back( 1 );
    }

    // FIXME This number seems to affect the value of the integration
    double scale = 1; 

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

    for ( unsigned int j = 0; j <= i; ++j ) {

        double r2 = rmesh[j];

    // Angular Integration
    double htot = 0; 
    for ( int ipart = 0; ipart < npart; ++ipart ) {
    
        double h_old = 0;
        double h_new = 0;
        double tol = 0.001;
        int ntheta = 64;
       
        // try with ntheta pts
        for( int m = 0; m < ntheta; ++m ) {

            // x = cos( theta )
            double x = ( xil.at(ipart) * ( ntheta - m - 1 ) 
                     + xiu.at(ipart) * m ) / ( ntheta - 1 );

            // integration weights
            double we = ( xiu.at(ipart) - xil.at(ipart) ) 
                      / ( ntheta - 1 );

            if( m == 0 ) we /= 2;
            if( m == ( ntheta - 1 ) ) we /= 2;

            double hs = volume_nl_form_PB1( r1, r2, x, R0, a0 );
            double hv = scale * gaussian_nl( r1, r2, x, beta );

                                      // angle in rad
            h_new += hs * hv * Pl.LegendreP0( l, std::acos( x ) ) * we; 

        } // end loop over angle

        // make sure integral converges
        while( std::abs( h_old - h_new ) > tol ) {

            h_old = h_new;

            double h = 0;
            ntheta *= 2;
            for( int m = 0; m < ntheta; ++m ) {

                // x = cos( theta )
                double x = ( xil[ipart] * ( ntheta - m - 1 ) + xiu[ipart] * m )
                         / ( ntheta - 1 );

                // integration weights
                double we = ( xiu[ipart] - xil[ipart] ) / ( ntheta - 1 );

                if( m == 0 ) we /= 2;
                if( m == ( ntheta - 1 ) ) we /= 2;

                double hs = volume_nl_form_PB1( r1, r2, x, R0, a0 );
                double hv = scale * gaussian_nl( r1, r2, x, beta );

                                          // angle in rad
                h += hs * hv * Pl.LegendreP0( l, std::acos( x ) ) * we; 

            } // end loop over angle
            h_new = h;
        } // end while loop
        htot += h_new;
    } // end ipart
        mtx( i, j ) = 2 * M_PI * htot / scale;
        if ( i != j ) mtx( j, i ) = 2 * M_PI * htot / scale;

    } // end loop over r2
    } // end loop over r1

    return mtx;
}

matrix_t 
q_volume_nl_PB1( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta, const LegendrePoly &Pl ) {

    // Prepare Angular Integration
    int npart;
    std::vector<double> xil;
    std::vector<double> xiu;

    if( l == 0 ) {
        npart = 1;
        xil.push_back( -1 );
        xiu.push_back( 1 );
    }
    else if ( l == 1 ) {
        npart = 2;
        xil.push_back( -1 );
        xiu.push_back( 0 );
        xil.push_back( 0 );
        xiu.push_back( 1 );
    }
    else if ( l == 2 ) {
        npart = 3;
        xil.push_back( -1 );
        xiu.push_back( -1 / std::sqrt( 3 ) );
        xil.push_back( -1 / std::sqrt( 3 ) );
        xiu.push_back( 1 / std::sqrt( 3 ) );
        xil.push_back( 1 / std::sqrt( 3 ) );
        xiu.push_back( 1 );

    }
    else {
        npart = 2;
        xil.push_back( -1 );
        xiu.push_back( 0 );
        xil.push_back( 0 );
        xiu.push_back( 1 );
    }

    // FIXME This number seems to affect the value of the integration
    double scale = 1; 
    double factor = std::pow( M_PI, 1.5 ) * std::pow( beta, 3 );

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

    for ( unsigned int j = 0; j <= i; ++j ) {

        double r2 = rmesh[j];

    // Angular Integration
    double htot = 0; 
    for ( int ipart = 0; ipart < npart; ++ipart ) {
    
        double h_old = 0;
        double h_new = 0;
        double tol = 0.001;
        int ntheta = 64;
       
        // try with ntheta pts
        for( int m = 0; m < ntheta; ++m ) {

            // x = cos( theta )
            double x = ( xil.at(ipart) * ( ntheta - m - 1 ) 
                     + xiu.at(ipart) * m ) / ( ntheta - 1 );

            // integration weights
            double we = ( xiu.at(ipart) - xil.at(ipart) ) 
                      / ( ntheta - 1 );

            if( m == 0 ) we /= 2;
            if( m == ( ntheta - 1 ) ) we /= 2;

            double q = std::sqrt( r1 * r1 + r2 * r2 + 2 * r1 * r2 * x );
            double hs = q_woods_saxon( q, R0, a0 );
            double hv = scale * gaussian_nl( r1, r2, x, beta );

                                      // angle in rad
            h_new += hs * hv * Pl.LegendreP0( l, std::acos( x ) ) * we; 

        } // end loop over angle

        // make sure integral converges
        while( std::abs( h_old - h_new ) > tol ) {

            h_old = h_new;

            double h = 0;
            ntheta *= 2;
            for( int m = 0; m < ntheta; ++m ) {

                // x = cos( theta )
                double x = ( xil[ipart] * ( ntheta - m - 1 ) + xiu[ipart] * m )
                         / ( ntheta - 1 );

                // integration weights
                double we = ( xiu[ipart] - xil[ipart] ) / ( ntheta - 1 );

                if( m == 0 ) we /= 2;
                if( m == ( ntheta - 1 ) ) we /= 2;

                double q = std::sqrt( r1 * r1 + r2 * r2 + 2 * r1 * r2 * x );
                double hs = q_woods_saxon( q, R0, a0 );
                double hv = scale * gaussian_nl( r1, r2, x, beta );

                                          // angle in rad
                h += hs * hv * Pl.LegendreP0( l, std::acos( x ) ) * we; 

            } // end loop over angle
            h_new = h;
        } // end while loop
        htot += h_new;
    } // end ipart
        mtx( i, j ) = 2 * M_PI * factor * htot / scale;
        if ( i != j ) mtx( j, i ) = 2 * M_PI * factor * htot / scale;

    } // end loop over r2
    } // end loop over r1

    return mtx;
}

// Calculates U * H, where U is the form given in the function
// volume_nl_form_PB2, and H is the Gaussian factor. 
matrix_t 
volume_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta ) {

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );

            double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                     * volume_nl_form_PB2( r1, r2, R0, a0 ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;
}
matrix_t 
q_volume_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
                 double a0, double beta ) {

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );
            double q = r1 + r2;

            double y = 4 * q_woods_saxon( q, R0, a0 ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;
}

matrix_t 
surface_nl_PB2( const std::vector<double> &rmesh, int l, double R0, 
                double a0, double beta ) {

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );

            double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                     * surface_nl_form_PB2( r1, r2, R0, a0 ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;

}

// Van Neck's form of nonlocality: sqrt[WS(r)WS(r^\prime)] * gaussian
matrix_t 
volume_nl_VN( const std::vector<double> &rmesh, int l, double R0, 
              double a0, double beta ) {
    

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );

            double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                     * volume_nl_form_VN( r1, r2, R0, a0 ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;
}

// nonlocal form for surface peaked potential ( in coordinate space )
matrix_t
surface_nl_VN( const std::vector<double> &rmesh, int l, double R0,
               double a0, double beta ) {

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );

            double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                     * surface_nl_form_VN( r1, r2, R0, a0 ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;
}

matrix_t 
q_form_bessel( const std::vector<double> &rmesh, int l, double R0, 
               double a0, double beta ) {
    

    // matrix to hold the nonlocal potential
    matrix_t mtx( rmesh.size(), rmesh.size() );

    // Do Angular Integration for each r1, r2
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        double r1 = rmesh[i];

        for ( unsigned int j = 0; j <= i; ++j ) {

            double r2 = rmesh[j];
            double x = 2 * r1 * r2 / std::pow( beta, 2 );

            double y = 4 / ( std::pow( beta, 3 ) * std::sqrt( M_PI ) ) 
                     * std::exp( - a0 * ( r1 - R0 ) ) 
                     * std::exp( - a0 * ( r2 - R0 ) ) 
                     * std::exp( - ( std::pow( ( r1 - r2 ) / beta, 2 ) ) )
                     * gsl_sf_bessel_il_scaled( l, x ); 

            mtx( i, j ) = y;
            if ( i != j ) mtx( j, i ) = y;

    } // end loop over r2
    } // end loop over r1


    return mtx;
}

// Step function
double step( double diff ) {
    
    if ( diff > 0. ) return 1.0;
    else return 0;
}

/*
    From "constants.h":

    m0 is the mass of the projectile nucleon
    e2 is the charge squared in MeV-fm
    kconstant is 2 / (\hbar **2)
*/


//relativistic conversion from Lab kinetic to total CM kinetic energy
double energyLab2Cm(double Elab, double A )
{
  if (Elab < 0.) return Elab*A/(1.+A);
 // center of mass velocity in units of c
 double vcm = std::sqrt(Elab*(Elab+2.*m0))/(Elab+(1.+A)*m0);
 //gamma factor for this velocity
 double gam = 1./std::sqrt(1.- std::pow(vcm,2));
 double Ecm = (gam-1.)*(1.+A)*m0 +
               gam*Elab*(A-1.)*m0/((A+1.)*m0+Elab);
 return Ecm;
}
//**********************************************************************
//relativistic conversion from CM kinetic to total Lab kinetic energy
double energyCm2Lab(double Ecm, double A )
{
  if (Ecm < 0.) return Ecm/A*(1.+A);

  //find momentum of projecile, also momentum of target
  double pc = std::sqrt(Ecm)*std::sqrt((Ecm+2.*m0)*(Ecm+2.*A*m0)*
              (Ecm+2.*(A+1.)*m0))/2./(Ecm+(A+1)*m0);
  //velocity of target in units of c
  double vtarget = pc/std::sqrt(std::pow(pc,2)+std::pow(A*m0,2));
  //gamma factor for this velocity
  double gam = 1./std::sqrt(1.- std::pow(vtarget,2));
  // tot energy of projectile (neutron or proton in com frame)
  double Eproj = std::sqrt(std::pow(m0,2)+std::pow(pc,2));
  double Elab = (Eproj + vtarget*pc)*gam;
  //this energy contains rest mass , so remove it 
  Elab -= m0;
  return Elab;
}

// Calculates the wave number given the total kinetic energy
// in the center of mass frame
double WaveNumber( double energyCM, double A, double EminRel ) { 

    // m0 is 931.5 MeV but is also used as the the nucleon mass.
    double mu_factor = A / ( A + 1 );
    double k2; //Wave number squared

    if ( energyCM <=  EminRel ) { 

        k2 = kconstant * mu_factor * energyCM;
    }
    else {
   
        double energyLab = energyCm2Lab( energyCM, A ); 

        k2 = kconstant * std::pow( A / ( 1. + A + energyCM/m0 ), 2 ) 
                       * energyLab * ( energyLab / 2. / m0 + 1. );

    }

    return std::sqrt( std::abs( k2 ) ); // Wave number in fm^-1

}

// Relativistic Correction
double Rel_Factor( double energyCM, double EminRel ) {

    if ( energyCM > EminRel ) 
        return 2 * ( energyCM / m0 + 1 ) / ( energyCM / m0 + 2 );
    else return 1;

}

//**************************************************************************
/* Subtracted Dispersion Integrals of the 'w' functions. Also included is 
   the dispersion integral of the asymmetric part of the volume imaginary 
   potential. 

   Parameters:
    A defines the magnitude
    B defines the rise away from the Fermi energy
    Ep gap width, \f$W=0\f$ for \f$|E-E_{fermi}| < E_{p} \f$
    Efermi is the Fermi Energy in MeV 
*/

// -------------------------------------------------------------------------
// Dispersion Integral of w( 2, ... ) from -Inf to Efermi
// -------------------------------------------------------------------------
double delta_w2_hole( double E, double Efermi, double A, 
                                               double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double EpB2 = std::pow( Ep, 2 ) + std::pow( B, 2 );

    if (Ex == 0.) return 0.;
    if (Ep == 0 ) return A*Ex/( std::pow(Ex,2)+ std::pow(B,2) )
                          *(B/2. + Ex * std::log( std::abs( Ex ) ) / pi );

    else if (Ex == -Ep) return -A*Ep/EpB2*( B/2.+Ep/pi*std::log( Ep/B ) );
    else { 
    
        double Eplus = Ex + Ep;
        return A/pi/EpB2 * ( std::pow( Ep, 2 ) 
                         * std::log( std::abs( Eplus/Ep ) )  
                         + ( Ex + 2. * Ep ) * Ex * std::pow( B, 2 ) 
                         / ( std::pow ( Eplus , 2 ) + std::pow( B, 2 ) ) 
                         * std::log( std::abs ( Eplus/B ) ) - Ex * B 
                         * ( Eplus * Ep - std::pow( B, 2 ) ) * pi/2.
                         / ( std::pow ( Eplus, 2) + std::pow( B, 2 ) ) );
    }
}

// -------------------------------------------------------------------------
// Dispersion Integral of w( 2, ... ) from Efermi to +Inf
// -------------------------------------------------------------------------
double delta_w2_particle(double E, double Efermi, double A, 
                                                  double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double EpB2 = std::pow( Ep, 2 ) + std::pow( B, 2 );

    if (Ex == 0.) return 0.;
    if (Ep == 0.) 
        return A*Ex/ ( std::pow(Ex,2) + std::pow(B,2) )
                *( B/2. - Ex * std::log( std::abs(Ex) )/pi );

    else if (Ex == Ep) return A*Ep/EpB2 *(B/2.+Ep/pi*log(Ep/B));
    else {
    
        double Eminus = Ex-Ep;
  
        return A/pi/EpB2 * (-std::pow( Ep, 2 ) 
                         * std::log( std::abs( Eminus/Ep ) )
                         - ( Ex - 2. * Ep ) * Ex * std::pow( B, 2 ) 
                         / ( std::pow( Eminus, 2 ) + std::pow( B, 2 ) ) 
                         * std::log( std::abs( Eminus/B ) ) + Ex * B 
                         * ( Eminus * Ep + std::pow(B,2) ) * pi/2.
                         / ( std::pow( Eminus, 2) + std::pow( B, 2 ) ) );
    }
}

// -------------------------------------------------------------------------
// Dispersion Integral of w( 4, ... ) from -Inf to Efermi
// -------------------------------------------------------------------------
double delta_w4_hole(double E, double Efermi, double A, 
                                              double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eplus = Ex + Ep;
	double SE1 = std::pow(Ep,4) + std::pow(B,4);
	double SE2 = std::pow(Eplus,4) + std::pow(B,4);

	double fact = std::sqrt(2.) * Ex / 4. 
                * ( B * std::pow( Ep * Eplus, 3 ) - std::pow( B, 7 ) ) 
                + ( std::pow( Eplus, 2 ) - std::pow( Ep, 2 ) ) / 4. 
                * ( std::pow( B, 6 ) - std::pow( B * Ep * Eplus, 2 ) ) 
                + ( std::pow( Eplus, 3) - std::pow( Ep, 3 ) ) *std::sqrt(2.)/4. 
                * ( std::pow( B, 3 ) * Ep * Eplus - std::pow( B, 5 ) );

	if (Eplus != 0.) fact -= std::pow( B * Eplus , 4 ) 
                           / pi * std::log( std::abs( Eplus/B ) );
                            
	if (Eplus != 0. && Ep != 0.) fact -= std::pow( Ep * Eplus, 4 ) 
                                       / pi* std::log( std::abs( Eplus/Ep ) );

	if (Ep != 0.) fact += std::pow( Ep * B, 4) / pi * std::log( Ep/B );

	return -A * fact / SE1 / SE2;
}

// -------------------------------------------------------------------------
// Dispersion Integral of w( 4, ... ) from Efermi to +Inf
// -------------------------------------------------------------------------
double delta_w4_particle(double E, double Efermi, double A, 
                                                  double B, double Ep ) {

    double pi = M_PI;
	double Ex = E - Efermi;
	double Eminus = Ex - Ep;
	double SE1 = std::pow(Ep,4) + std::pow(B,4);
	double SE2 = std::pow(Eminus,4) + std::pow(B,4);
	 
	double fact = ( std::pow(B,7) + B * std::pow( Ep * Eminus, 3 ) ) 
                * Ex * std::sqrt(2.) / 4.  
                + ( std::pow( Eminus, 2 ) - std::pow( Ep, 2 ) )
                * (std::pow( B, 6 ) - std::pow( B * Ep * Eminus, 2 ) ) / 4. 
                + (std::pow( Eminus, 3 ) + std::pow( Ep, 3 ) ) * std::sqrt(2.)
                * (std::pow( B, 5 ) + pow( B, 3 ) * Ep * Eminus ) / 4.;

	if (Eminus != 0. && Ep!= 0.) fact -= std::pow( Ep * Eminus, 4 ) / pi
                                       * std::log( std::abs( Eminus/Ep ) );

	if (Eminus !=0) fact -= std::pow( B * Eminus, 4 ) / pi 
                          * std::log(std::abs( Eminus/B ) ) ;

	if (Ep != 0.) fact += std::pow( Ep * B, 4 ) / pi * std::log( Ep/B );

	return A * fact / SE1 / SE2; 
}

// -------------------------------------------------------------------------
// Dispersion Integral of asymmetric correction from -Inf to Efermi 
// -------------------------------------------------------------------------
double delta_asymmetric_hole( double E, double Efermi, double A, double Ea ) {


    double pi = M_PI;
    double Ex = E - Efermi;

    if (Ex == 0.) return 0.;
    
    if ( Ea == 0 ) return -A * Ex / ( std::pow( Ex, 2 ) + std::pow( Ea, 2 ) )
                          * ( Ea / 2. + Ex * std::log( std::abs( Ex ) ) / pi );
    else if ( Ex == -Ea ) return A / 4.;
    else { 
        
	    double Eplus = Ex + Ea;
        return -A /pi/2. * ( std::log( std::abs( Eplus/Ea ) ) 
               + ( Ex + 2. * Ea ) * Ex 
               / ( std::pow( Eplus, 2 ) + std::pow( Ea, 2 ) ) 
               * std::log( std::abs( Eplus/Ea ) )
               - std::pow( Ex, 2 ) * pi / 2. 
               / ( std::pow( Eplus, 2 ) + std::pow( Ea, 2 ) ) ); 
    }
}

// -------------------------------------------------------------------------
// Dispersion Integral of asymmetric correction from Efermi to +Inf
// -------------------------------------------------------------------------
double delta_asymmetric_particle(double E, double Efermi, 
                                           double alpha, double Ea) {

    double pi = M_PI;
    double El = Efermi + Ea;
    double term;

    if (E == El) {
    
        term = std::sqrt( El ) / 2. 
             * ( std::log( 16. ) + 3. * std::log( std::abs( El/Ea ) ) );

        term -= std::pow( El, 3./2. ) / 2. 
              / Efermi * std::log( std::abs( El/Ea ) );

        term += 2.* std::sqrt( std::abs( Efermi ) ) 
              * ( pi/2. - std::atan( std::sqrt( std::abs( El/Efermi ) ) ) );

        return term*alpha/pi;
    }  

    if (E < 0.) term = -2. * std::sqrt( std::abs( E ) ) 
                     * ( pi/2. - std::atan( std::sqrt( std::abs( El/E ) ) ) );

    else term = std::sqrt( E ) * std::log( std::abs( ( std::sqrt( El ) 
              + std::sqrt( E ) ) / ( std::sqrt( El ) - std::sqrt( E ) ) ) );

    term += 2. * std::sqrt( std::abs( Efermi ) ) 
          * ( pi/2.- std::atan( std::sqrt( El / std::abs( Efermi ) ) ) );

    if (E == 0.) term += std::pow( El, 3./2. ) / 2. 
                       *( 1./El - std::log( std::abs( El/Ea ) ) / Efermi );

    else term += std::pow( El, 3./2. ) / 2. / E / Efermi 
               * ( Efermi * std::log( std::abs( El/(El-E) ) ) 
               - E * std::log( std::abs( El/Ea ) ) );

    term += 3./2. * std::sqrt( El ) * std::log( std::abs( (El-E)/Ea ) );

    return term * alpha / pi;

}

//***********************************************************************
  //principle value integral for m =2,C==0 from analytic expression of
  //VanderKam, J. Phys. G. 26 (2000) 1787
double delta_w2( double E, double Efermi, double A, double B, double Ep ) {

    double Ex = E - Efermi;
    double pi = M_PI;

    if ( Ex == 0. ) return 0.;
    else if ( Ep == 0 ) return 
        A * B * Ex / ( std::pow( Ex, 2 ) + std::pow( B, 2 ) );

    else if ( std::abs( Ex ) == Ep ) return 
        A / pi / ( std::pow( Ex + Ep, 2 ) + std::pow( B, 2 ) ) 
          / ( std::pow( Ex - Ep, 2) + std::pow( B, 2) ) 
          * ( pi * B * Ex * ( Ex * Ex - Ep * Ep + B * B ) 
                 + 2. * std::pow( B * Ep, 2 ) 
                 * std::log( 4. * std::pow( Ep / B, 2 ) ) );

    else return A / pi / ( std::pow( Ex + Ep, 2 ) + std::pow( B, 2 ) ) 
        / ( std::pow( Ex - Ep, 2 ) + std::pow( B, 2 ) ) 
        * ( pi * B * Ex * ( Ex * Ex - Ep * Ep + B * B ) 
               + ( std::pow( Ex * Ex - Ep * Ep, 2 ) 
                   + B * B * ( Ex * Ex + Ep * Ep ) ) 
               * std::log( std::abs( ( Ex + Ep ) / ( Ex - Ep ) ) ) 
               + 2. * B * B * Ex * Ep 
               * std::log( std::abs( ( Ex * Ex - Ep * Ep ) / B / B ) ) );

}
//***********************************************************************
  //principle value integral for m =4,C==0 from analytic expression of
  //VanderKam, J. Phys. G. 26 (2000) 1787
double delta_w4( double E, double Efermi, double A, double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    if (Ex == 0.) return 0.;

    double fact = pi * B * B * Ex * Ep 
        * ( std::pow( Ex * Ex - Ep * Ep, 2 ) - std::pow( B, 4 ) ) 
        + pi * B * Ex / std::sqrt( 2. ) 
        * ( std::pow( Ex * Ex - Ep * Ep, 3 ) + B * B 
            * ( std::pow( B, 4 ) + 3. * std::pow( B * Ep, 2 ) 
                + std::pow( B * Ex, 2 ) + std::pow( Ex, 4 ) 
                + 2. * std::pow( Ex * Ep, 2 ) - 3. * std::pow( Ep, 4 ) ) );

    if ( Ex == Ep ) fact += 8. * std::pow( B * Ep, 4 ) 
        * std::log( 4. * std::pow( Ep / B, 2 ) );
    else if ( Ex == -Ep ) fact -= 8. * std::pow( B * Ep, 4 ) 
        * std::log( 4. * std::pow( Ep / B, 2 ) );
    else if ( Ep != 0. ) fact += 4. * std::pow( B, 4 ) 
        * Ex * Ep * ( Ex * Ex + Ep * Ep ) 
        * std::log( std::abs( ( Ex * Ex - Ep * Ep ) / B / B ) ) 
        + ( std::pow( Ex * Ex - Ep * Ep, 4 ) + std::pow( B, 4 ) 
            * ( std::pow( Ep, 4 ) + 6. * std::pow( Ep * Ex, 2 ) 
                + std::pow( Ex, 4 ) ) ) 
        * std::log( std::abs( ( Ex + Ep ) / ( Ex - Ep ) ) );

    return  A / pi / ( std::pow( Ex + Ep, 4 ) + std::pow( B, 4 ) ) 
              / ( std::pow( Ex - Ep, 4 ) + std::pow( B, 4 ) ) * fact;

  /*
  return A/pi/(pow(Ex+Ep,4)+pow(B,4))/(pow(Ex-Ep,4)+pow(B,4))*
    (pi*B*B*Ex*Ep*(pow(Ex*Ex-Ep*Ep,2)-pow(B,4))    +
       pi*B*Ex/sqrt(2.)*(pow(Ex*Ex-Ep*Ep,3)+B*B*(pow(B,4)+3.*pow(B*Ep,2)
       + pow(B*Ex,2)+pow(Ex,4)+2.*pow(Ex*Ep,2)-3.*pow(Ep,4)))   +
       4.*pow(B,4)*Ex*Ep*(Ex*Ex+Ep*Ep)*log(abs((Ex*Ex-Ep*Ep)/B/B)) +
       (pow(Ex*Ex-Ep*Ep,4)+pow(B,4)*(pow(Ep,4)+6.*pow(Ep*Ex,2)+pow(Ex,4)))*
       log(abs((Ex+Ep)/(Ex-Ep))));
  */
}

//**************************************************************************
/* Derivatives of the subtracted dispersion relations
*/

// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of w2_hole
// -------------------------------------------------------------------------
double der_delta_w2_hole( double E, double Efermi, double A, 
                                                   double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eplus = Ex + Ep;
    double EB = std::pow( Eplus, 2 ) + std::pow( B, 2 );

    if (Ex == -Ep) return A/2./B;
    else if (Ep == 0.) 
        return A * ( std::pow( B, 2 ) - std::pow( Ex, 2 ) ) / std::pow( EB, 2 ) 
                 * ( B/2. + Ex/pi * std::log( std::abs( Ex ) ) ) 
                 + A * Ex/pi/EB * ( 1. + std::log( std::abs( Ex ) ) );

    return A/2./pi / std::pow( EB, 2 ) * ( 2. * Eplus * EB 
                   + B * ( std::pow( B, 2 ) - std::pow( Eplus, 2 ) ) * pi 
			       + 4. * std::pow( B, 2 ) * Eplus 
                   * std::log( std::abs( Eplus/B ) ) );
}

// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of w2_particle
// -------------------------------------------------------------------------
double der_delta_w2_particle( double E, double Efermi, double A,
                                                       double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eminus = Ex - Ep;
    double EB = std::pow( Eminus, 2 ) + std::pow( B, 2 );

    if (Ex == Ep) return A/2./B;
    else if (Ep == 0.) 
        return A * ( std::pow( B, 2 ) - std::pow( Ex, 2 ) ) / std::pow( EB, 2 )
                 * ( B/2. - Ex/pi * std::log( std::abs( Ex ) ) ) 
                 - A * Ex/pi / EB * ( 1. + std::log( std::abs( Ex ) ) );

    return A/2./pi / std::pow( EB, 2 ) * ( -2. * Eminus * EB 
                   + B * ( std::pow( B, 2 ) - std::pow( Eminus, 2 ) ) * pi 
			       - 4. * std::pow( B, 2 ) * Eminus 
                   * std::log( std::abs( Eminus/B ) ) );
}

// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of w4_hole
// -------------------------------------------------------------------------
double der_delta_w4_hole( double E, double Efermi, double A, 
                                                   double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eplus = Ex + Ep;
    double SE1 = std::pow( Ep, 4 ) + std::pow( B, 4 );
    double SE2 = std::pow( Eplus, 4 ) + std::pow( B, 4 );

    if (Ex == -Ep) return A/2. / std::sqrt( 2. ) / B;

    else if (Ep == 0.) 
        return A * ( B * ( B - Ex ) * ( B + Ex ) 
                 * ( -2. * B * Ex * ( std::pow( B, 2 ) + std::pow( Ex,2 ) ) 
                 + std::sqrt( 2. ) * ( std::pow( B, 4 ) 
                 + 4. * std::pow( B,2 ) * std::pow( Ex, 2 ) 
                 + std::pow( Ex, 4 ) ) ) * pi  
                 + 4. * std::pow( Ex, 3 ) * ( SE2 + 4. * std::pow( B, 4 ) 
                 * std::log( std::abs( Ex/B ) ) ) ) 
                 / ( 4. * std::pow( std::pow( B, 4 ) 
                 + std::pow( Ex, 4 ), 2 ) * pi ); 

    else return A * ( SE1 * ( 4. * std::pow( Eplus,3 ) * SE2 
                  + B * ( B - Eplus ) * ( B + Eplus ) 
                  * ( -2. * B * Eplus * (std::pow( B, 2 ) 
                  + std::pow( Eplus, 2 ) ) + std::sqrt( 2. ) 
                  * ( std::pow( B, 4 ) + 4. * std::pow( B, 2 ) 
                  * std::pow( Eplus, 2 ) + std::pow( Eplus, 4 ) ) ) * pi ) 
                  + 16. * std::pow( B, 4 ) * std::pow( Eplus, 3 ) 
                  * ( std::pow( B, 4 ) * std::log( std::abs( Eplus ) / B ) 
                  + std::pow( Ep, 4 ) * ( std::log( Ep/B ) 
                  + std::log( std::abs( Eplus )/Ep ) ) ) ) 
                  / ( 4. * SE1 * std::pow( SE2, 2) * pi );

}

// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of w4_particle
// -------------------------------------------------------------------------
double der_delta_w4_particle( double E, double Efermi, double A,
                                                       double B, double Ep ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eminus = Ex - Ep;
    double SE2 = std::pow( Eminus, 4 ) + std::pow( B, 4 );
    double SE1 = std::pow( Ep, 4 ) + std::pow( B, 4 );

    if (Ex == Ep) return A/2./ std::sqrt( 2. ) / B;

    else if (Ep == 0.) 
        return -A * ( 4. * std::pow( B, 4 ) * std::pow( Ex, 3 ) 
                  + 4. * std::pow( Ex, 7 ) - std::sqrt( 2. ) 
                  * std::pow( B, 7 ) * pi - 2. * std::pow( B, 6 ) * Ex * pi 
                  - 3. * std::sqrt( 2. ) * std::pow( B, 5 ) * std::pow( Ex, 2) 
                  * pi + 3. * std::sqrt( 2. ) * std::pow( B, 3 ) 
                  * std::pow( Ex, 4 ) * pi + 2. * std::pow( B, 2 ) 
                  * std::pow( Ex, 5 ) * pi + std::sqrt( 2. ) * B 
                  * std::pow( Ex, 6 ) * pi + 16. * std::pow( B, 4 ) 
                  * std::pow( Ex, 3 ) * std::log( std::abs( Ex/B ) ) ) 
                  / ( 4. * std::pow( std::pow( B, 4 ) 
                  + std::pow( Ex, 4 ), 2 ) * pi );
  
  else return A * ( B * SE1 * ( std::sqrt(2.) * ( std::pow( B, 4 )  
                + 4. * std::pow( B, 2 ) * std::pow(Eminus,2) 
                + std::pow( Eminus, 4 ) ) + 2. * B * ( std::pow( B, 2 ) 
                + std::pow( Eminus, 2 ) ) * Eminus ) 
                * ( std::pow( B, 2 ) - std::pow( Eminus, 2 ) ) * pi 
                - 4. * std::pow( Eminus, 3 ) * ( SE1 * SE2  
                + 4. * std::pow( B, 4 ) * ( std::pow( B, 4 ) 
                * std::log( std::abs( Eminus/B ) ) 
                + std::pow( Ep, 4 ) * ( std::log( Ep/B ) 
                + std::log( std::abs( Eminus/Ep ) ) ) ) ) ) 
                / ( 4. * ( std::pow( B, 4 ) + std::pow( Ep, 4 ) ) 
                * std::pow( SE2, 2 ) * pi);

}

// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of asymmetic_hole
// -------------------------------------------------------------------------
double der_delta_asymmetric_hole( double E, double Efermi, double A, double Ea ) {

    double pi = M_PI;
    double Ex = E - Efermi;
    double Eplus = Ex + Ea;
    double EB = std::pow( Eplus, 2 ) + std::pow( Ea, 2 );

    if (Ex == -Ea) return -A/2./Ea;

    return -A/2./pi / std::pow( EB, 2 ) * ( 2. * Eplus * EB 
                    + Ea * ( std::pow( Ea, 2 ) - std::pow( Eplus, 2) ) * pi 
			        + 4. * std::pow( Ea, 2 ) * Eplus 
                    * std::log( std::abs( Eplus/Ea ) ) );
}
// -------------------------------------------------------------------------
// Derivative of Dispersion Integral of asymmetric_particle
// -------------------------------------------------------------------------
double der_delta_asymmetric_particle( double E, double Efermi, 
                                                double alpha, double Ea ) {

    double pi = M_PI;
    double El = Efermi + Ea;

    if (E == El) return alpha * ( 1. + std::log(4.) ) / 2. / pi / std::sqrt(El);
    else if ( std::abs(E) < 0.01 )  
        return 3. * alpha / 4. / std::sqrt(El) / pi
                  - 5./2. * E / std::pow( El, 3./2. ) / pi * alpha;

    else if (E > 0.)   
        return alpha/2./pi / std::pow( E, 2 ) * (std::pow( E, 3./2. ) 
                           * std::log( ( std::sqrt(El) + std::sqrt(E) )  
                           / std::abs( std::sqrt(El) - std::sqrt(E) ) ) 
                           + std::sqrt(El) * ( E  
                           - El * std::log( std::abs( El/( El - E ) ) ) ) );

    else return alpha/2./pi * ( 2. * std::sqrt( std::abs(El) ) 
                     / ( std::abs(El) - E ) + ( pi 
                     - 2. * std::atan( std::sqrt( std::abs( El/E ) ) ) ) 
                     / std::sqrt( std::abs(E) ) + std::sqrt(El) 
                     * ( ( El - 3. * E ) * E + El * ( E - El ) 
                     * std::log( El/( El - E ) ) ) 
                     / std::pow( E, 2 ) / ( El - E ) );

}

