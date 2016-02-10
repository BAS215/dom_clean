// functions related to calculating the charge density and matter density

#include "density.h"

// d_mtx_lj is the density matrix in coordinate space 
// with orbital angular momentum l and total angular momentum j. 
std::vector<double>
point_density_lj( const std::vector<double> &rmesh, const matrix_t &d_mtx_lj,
                  double xj ) {

    std::vector<double> point_dist;

    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        point_dist.push_back( ( 2 * xj + 1 ) * d_mtx_lj( i, i ) 
                            / ( 4 * M_PI ) );
    }

    return point_dist;
}

// Contribution to the point density distribution from the quasiholes
// The eigenfunction QPF should be normalized to 1
std::vector<double>
point_density_qp_nlj( const std::vector<double> &QPF, double xj ) {

    std::vector<double> point_dist_qp;

    for ( unsigned int i = 0; i < QPF.size(); ++i ) {

        point_dist_qp.push_back( ( 2 * xj + 1 ) * QPF[i] * QPF[i] 
                               / ( 4 * M_PI ) );
    }

    return point_dist_qp;
}

// relativistic spin-orbit correction
std::vector<double>
so_correction_nlj( const std::vector<double> &rmesh, 
                   const std::vector<double> &chd_nlj, 
                   int l, double xj, double tz, double occupation ) {

    double lp;
    if ( xj > l ) lp = static_cast<double>( l ); 
    else lp = static_cast<double>( - l - 1 );

    double mag_moment;
    if ( tz > 0 ) mag_moment = 2.29; // actually mu_p - 0.5
    else mag_moment = -1.91; // mu_n;

    std::vector<double> so_correction;
    for( unsigned int i = 1; i < rmesh.size() - 1; ++i ) {
        
        // derivative of point density distribution times r
        double rdelt = std::abs( rmesh[i] - rmesh[i-1] );
        double der = der_3pt( rdelt, rmesh[i-1] * chd_nlj[i-1],
                              rmesh[i+1] * chd_nlj[i+1] );
                              
        double term = 0.5 * occupation * lp * mag_moment * der
                    / std::pow( rmesh[i], 2 );

        so_correction.push_back( term );
    }

    return so_correction;
}

// point distribution folded with nucleon charge distribution
std::vector<double>
folded_ch_density( const std::vector<double> &rmesh, 
                   const std::vector<double> &rweights,
                   const std::vector<double> &point_dist, 
                   double tz, double A ) {
    
    // Experimental density parameterized as sum of gaussians: 
    // G( r ) = \Sum_{i}(theta[i] * exp(-r^2/r[i]^2) / (pi * r[i]^2)^(3/2) )
    // The numbers are from J. Phys. G: Nucl. Phys. 5, 1655 ( 1979 )
    std::vector< double > theta_g; // theta parameters
    std::vector< double > r_g2; // radial parameters (squared)

    // 0.21 fm is hbar_over_mc 
    double cm_correction = 0.5 * 0.21 * 0.21  - std::pow( A, - 2 / 3.0 );
    if( tz > 0 ) { 
        
        theta_g.push_back( 0.506 );
        theta_g.push_back( 0.328 );
        theta_g.push_back( 0.166 );
        r_g2.push_back( 0.432 + cm_correction );
        r_g2.push_back( 0.139 + cm_correction );
        r_g2.push_back( 1.526 + cm_correction );
    }
    else {

        theta_g.push_back( 1 );
        theta_g.push_back( -1 );
        r_g2.push_back( 0.469 + cm_correction );
        r_g2.push_back( 0.546 + cm_correction );

    }

    //loop over gaussians
    std::vector<double> chdf;
    chdf.assign( rmesh.size(), 0 );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        double r = rmesh[i];
        for( unsigned int k = 0; k < r_g2.size(); ++k ) {

            double r_g = std::sqrt( r_g2[k] );

            // integral
            for( unsigned int j = 0; j < rmesh.size(); ++j ) {

                double rp = rmesh[j];
                chdf[i] += rweights[j] * ( theta_g[k] / r_g ) * rp / r 
                         * ( std::exp( - std::pow( ( r - rp ) / r_g, 2 ) )
                           - std::exp( - std::pow( ( r + rp ) / r_g, 2 ) ) )
                         * point_dist[j] / std::sqrt( M_PI );
            } // end integral
        } // end sum

    }

    return chdf;
   // return point_dist;

}

// point distribution folded with finite size of nucleon 
// finite size of nucleon for both neutrons and protons is
// approximated with the proton charge distribution
std::vector<double>
matter_distribution( const std::vector<double> &rmesh, 
                     const std::vector<double> &rweights,
                     const std::vector<double> &point_dist, 
                     double A ) {
    
    // Experimental density parameterized as sum of gaussians: 
    // G( r ) = \Sum_{i}(theta[i] * exp(-r^2/r[i]^2) / (pi * r[i]^2)^(3/2) )
    // The numbers are from J. Phys. G: Nucl. Phys. 5, 1655 ( 1979 )
    std::vector< double > theta_g; // theta parameters
    std::vector< double > r_g2; // radial parameters (squared)

    theta_g.push_back( 0.506 );
    theta_g.push_back( 0.328 );
    theta_g.push_back( 0.166 );
    r_g2.push_back( 0.432 );
    r_g2.push_back( 0.139 );
    r_g2.push_back( 1.526 );

    //loop over gaussians
    std::vector<double> chdf;
    chdf.assign( rmesh.size(), 0 );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        double r = rmesh[i];
        for( unsigned int k = 0; k < r_g2.size(); ++k ) {

            double r_g = std::sqrt( r_g2[k] );

            // integral
            for( unsigned int j = 0; j < rmesh.size(); ++j ) {

                double rp = rmesh[j];
                chdf[i] += rweights[j] * ( theta_g[k] / r_g ) * rp / r 
                         * ( std::exp( - std::pow( ( r - rp ) / r_g, 2 ) )
                           - std::exp( - std::pow( ( r + rp ) / r_g, 2 ) ) )
                         * point_dist[j] / std::sqrt( M_PI );
            } // end integral
        } // end sum

    }

    return chdf;

}
