//
#include "read_data.h"
#include "density.h"

int main ( ) {

    Data_set D = read_general_set( "exp_chd_ca40_sog.out" );

    std::vector< double > fguess;
    std::vector< double > rweights;
    double deltaR = D.x_values.at( 1 ) - D.x_values.at( 0 );
    for ( unsigned int i = 0; i < D.x_values.size(); ++i ) {

        fguess.push_back( D.y_values[i] );
        rweights.push_back( deltaR );
    }

    // First iteration
    double tz = 0.5;
    double A = 40;

    // Convolution with finite size of proton charge density (tz = 0.5)
    std::vector< double > pconvoluted = 
        folded_ch_density( D.x_values, rweights, fguess, 0.5, A );

    // Convolution with finite size of neutron charge density (tz = -0.5)
    std::vector< double > nconvoluted = 
        folded_ch_density( D.x_values, rweights, fguess, -0.5, A );

    // First ratio
    std::vector< double > pcnt_ratio;
    std::vector< double > fconvoluted; 
    for ( unsigned int i = 0; i < D.x_values.size(); ++i ) {

        // add neutron and proton convolutions together
        fconvoluted.push_back( pconvoluted[i] + nconvoluted[i] );

        // calculate ratio
        double r = D.y_values[i] / fconvoluted[i];

        // calculate percent difference
        pcnt_ratio.push_back( std::abs( r - 1.00 ) );

        // calculate new guess
        fguess[i] = fguess[i] * r;

    }

    double max_pcnt = 0;
    for ( unsigned int i = 0; i < pcnt_ratio.size(); ++i ) {

        max_pcnt = std::max( pcnt_ratio[i], max_pcnt ); 
    }
/*
    // Write out first iteration to file
    std::ofstream tfile( "unfolded_first.out" );
    fconvoluted.clear();
    fconvoluted = folded_ch_density( D.x_values, rweights, fguess, tz, A );
    for ( unsigned int i = 0; i < fguess.size(); ++i ) {

        tfile << D.x_values[i] << " " << fguess[i] 
                               << " " << fconvoluted[i] 
                               << " " << D.y_values[i] << std::endl;
    }

    tfile.close();
    tfile.clear();
*/

    double tolerance = 0.01;
    int max_iterations = 100;
    int counter = 0;
    while ( max_pcnt > tolerance ) {

        fconvoluted.clear();
        pcnt_ratio.clear();

        // Convolution with finite size of proton charge density (tz = 0.5)
        pconvoluted.clear();
        pconvoluted = folded_ch_density( D.x_values, rweights, fguess, 0.5, A );

        // Convolution with finite size of neutron charge density (tz = -0.5)
        nconvoluted.clear();
        nconvoluted = folded_ch_density( D.x_values, rweights, fguess, -0.5, A);

        for ( unsigned int i = 0; i < D.x_values.size(); ++i ) {

            // add neutron and proton convolutions together
            fconvoluted.push_back( pconvoluted[i] + nconvoluted[i] );

            double r = D.y_values[i] / fconvoluted[i];

            pcnt_ratio.push_back( std::abs( r - 1.00 ) );

            fguess[i] = fguess[i] * r;

        }

        max_pcnt = 0;
        for ( unsigned int i = 0; i < pcnt_ratio.size(); ++i ) {

            max_pcnt = std::max( pcnt_ratio[i], max_pcnt ); 
        }

        counter++;

        if ( counter >= max_iterations ) {

            std::cout << "maximum number of iterations reached." << std::endl;
            break;
        }
    }

    std::cout << "max_pcnt = " << max_pcnt * 100 << std::endl;
    std::cout << "Number of Iterations = " << counter << std::endl;

    // Write out results to file
    std::ofstream file( "unfolded.out" );

    fconvoluted.clear();
    pconvoluted.clear();
    nconvoluted.clear();

    pconvoluted = folded_ch_density( D.x_values, rweights, fguess, 0.5, A );
    nconvoluted = folded_ch_density( D.x_values, rweights, fguess, -0.5, A);

    double norm = 0;
    double chd_r2 = 0;
    for ( unsigned int i = 0; i < fguess.size(); ++i ) {

        fconvoluted.push_back( pconvoluted[i] + nconvoluted[i] );

        file << D.x_values[i] << " " << fguess[i] 
                              << " " << fconvoluted[i] 
                              << " " << D.y_values[i] 
                              << " " << pcnt_ratio[i] << std::endl;

        norm += 4 * M_PI * rweights[i] * D.x_values[i] * D.x_values[i]
              * fconvoluted[i];

        chd_r2 += 4 * M_PI * rweights[i] * std::pow( D.x_values[i], 4 ) 
                * fconvoluted[i];
    }

    double rms = std::sqrt( chd_r2 / norm );

    std::cout << "Rrms = " << rms << std::endl;
    std::cout << "Norm = " << norm << std::endl;

    file.close();
    file.clear();

    return 0;
}
