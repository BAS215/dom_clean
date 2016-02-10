//
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include "legendre.h"
#include "boundRspace.h"
#include "meshes.h"
#include "numerical.h"
#include "io.h"
#include "Gauss_q.h"


int main( ) {

    std::string input_dir = "Input/";

    // Read Nucleus File
    std::string nucleus_string;
    std::ifstream nucleus_file( "domgen.config" );

    nucleus_file >> nucleus_string;

    nucleus_file.close();
    nucleus_file.clear();

    // Read Configuration file
    std::string config_filename = nucleus_string + ".config";
    std::ifstream config_file( config_filename.c_str() );

    std::string parameters_string;
    int fit_ph;
    double rmax;
    int rpts;
    int lmax;

    config_file >> parameters_string >> fit_ph;
    config_file >> rmax >> rpts >> lmax;

    int num_lj;
    config_file >> num_lj;

    config_file.close();
    config_file.clear();

    // Read in normalization
    double calc_Z = 20;
    double calc_N = 20;
    std::cout << "Enter calculated Z: " << std::endl;
    std::cin >> calc_Z; 

    std::cout << "Enter calculated N: " << std::endl;
    std::cin >> calc_N;

    std::string output_dir = "Output_" + parameters_string + "/Spectral_functions_kE/";
    std::string parameters_filename = parameters_string + ".inp";

    std::cout << "rmax = " << rmax << std::endl;
    std::cout << "rpts = " << rpts << std::endl;
    std::cout << "lmax = " << lmax << std::endl;

    std::string n_string = "n" + nucleus_string;
    std::string p_string = "p" + nucleus_string;

    std::string n_filename = input_dir + n_string + ".inp";
    std::string p_filename = input_dir + p_string + ".inp";

    // Create Nuclear Parameter Objects
    NuclearParameters Nu_n = read_nucleus_parameters( n_filename );
    NuclearParameters Nu_p = read_nucleus_parameters( p_filename );

    std::vector< NuclearParameters > Nu_vec;
    Nu_vec.push_back( Nu_n );
    Nu_vec.push_back( Nu_p );


    // Read in DOM parameters
    std::ifstream pfile( parameters_filename.c_str() );
    if ( pfile.is_open() !=1 ) {
        std::cout << "could not open file " << parameters_filename << std::endl;
        std::abort();
    }
    pfile.close();
    pfile.clear();

    // Create radial grid
    std::vector<double> rmesh;
    std::vector<double> rweights;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        rweights.push_back( rdelt );
        
    }

    // Create momentum space grid for k-slice
    std::vector<double> kmesh1;
    double kmax1 = 6.0;
    int kpts1 = 100;
    double deltak1 = kmax1 / kpts1;

    for ( int i = 0; i < kpts1; ++i ) {

        kmesh1.push_back( ( i + 0.5 ) * deltak1 );
    }

    // Create momentum space grid for E-slice
    std::vector<double> kmesh2; // wave number in units of fm^{-1}
    std::vector<double> pmesh; // momentum in units of MeV / c
    double pmin = 170; 
    double pmax = 650;
    double deltap = 40;
    int ppts = static_cast<int> ( ( pmax - pmin ) / deltap ) + 1;

    for ( int i = 0; i < ppts; ++i ) {

        double p = pmin + i * deltap;
        pmesh.push_back( p );

        double k = p / hbarc;
        kmesh2.push_back( k );
    }

    // Prepare stuff for output files
    std::vector< std::string > np_strings;
    np_strings.push_back( n_string );
    np_strings.push_back( p_string );

    // create L and J strings for output files
    std::string L_array[8] = { "s", "p", "d", "f", "g", "h", "i", "j" };
    std::string J_array[8] = { "1", "3", "5", "7", "9", "11", "13", "15" };

    //Potential specifications
    int type = 1; // 1 is P.B. form (average), 0 is V.N. form
    int mvolume = 4;
    int AsyVolume = 1;


    /* CALCULATIONS */

    // Loop over protons and neutrons
    double Zp; // number of protons in projectile
    for ( unsigned int nu = 0; nu < Nu_vec.size(); ++nu ) {
        
        double tz = nu - 0.5; // +0.5 for protons, -0.5 for neutrons
        double calc_norm;
        if ( tz < 0 ) {

            Zp = 0;
            calc_norm = calc_N;
        }
        else { 

            Zp = 1;
            calc_norm = calc_Z;
        }

        // Nucleus Object
        const NuclearParameters &Nu = Nu_vec[nu];

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

        // Construct Potential Object
        pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
        pot *U1 = &U;

        // Construct Object for calculating bound-state properties
        boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, lmax, Nu.Z, Zp, Nu.A, U1 );

        // Create Energy Grid for k-slice
        double Emin1 = -225;
        double Emax1 = -25;
        double deltaE1 = 25;
        int epts1 = static_cast<int>( ( Emax1 - Emin1 ) / deltaE1 ) + 1;
        std::vector<double> emesh1;
        for ( int i = 0; i < epts1; ++i ) {

            emesh1.push_back( Emin1 + i * deltaE1 );
        }

        // Create Energy grid for E-slice
        double Emin2 = -300;
        double Emax2 = -25;
        double deltaE2 = 2;
        int epts2 = static_cast<int>( ( Emax2 - Emin2 ) / deltaE2 ) + 1;
        std::vector<double> emesh2;
        for ( int i = 0; i < epts2; ++i ) {

            emesh2.push_back( Emin2 + i * deltaE2 );
        }

        matrix_t S_of_kE_mtx1( kmesh1.size(), emesh1.size() );
        matrix_t S_of_kE_mtx2( kmesh2.size(), emesh2.size() );

        // intialize matrices to zero
        for ( unsigned int i = 0; i < kmesh1.size(); ++i ) {
            for ( unsigned int j = 0; j < emesh1.size(); ++j ) {

                S_of_kE_mtx1( i, j ) = 0;
            }
        }

        for ( unsigned int i = 0; i < kmesh2.size(); ++i ) {
            for ( unsigned int j = 0; j < emesh2.size(); ++j ) {

                S_of_kE_mtx2( i, j ) = 0;
            }
        }

        std::vector< matrix_t > S_of_kE_mtx2_lj_vec;

        // Loop over lj channels
        for ( int L = 0; L < lmax + 1; ++L ) {

            for( int up = -1; up < 2; up+=2 ) {

                double xj = L + up / 2.0;
                int j_index = ( up - 1 ) / 2;
                if ( xj < 0 ) continue;

                std::string j_string;
                if ( L == 0 ) j_string = J_array[ L ];
                else j_string = J_array[ L + j_index ];

                std::string lj_string = L_array[L] + j_string + "2";

                // Create Bessel Function matrix in k and r
                matrix_t bess_mtx1( kmesh1.size(), rmesh.size() );
                for( unsigned int nk = 0; nk < kmesh1.size(); ++nk ) {
                for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

                    double rho = kmesh1[nk] * rmesh[nr];

                    bess_mtx1( nk, nr ) = gsl_sf_bessel_jl( L, rho );
                }
                }

                // Create Bessel Function matrix in k and r
                matrix_t bess_mtx2( kmesh2.size(), rmesh.size() );
                for( unsigned int nk = 0; nk < kmesh2.size(); ++nk ) {
                for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

                    double rho = kmesh2[nk] * rmesh[nr];

                    bess_mtx2( nk, nr ) = gsl_sf_bessel_jl( L, rho );
                }
                }

                // Calculate S( k; E ) for k-slice 
                for ( unsigned int m = 0; m < emesh1.size(); ++m ) {
            
                    double E = emesh1[m];

/*
                    std::ostringstream e_strm;
                    e_strm << "m" << std::abs( E );
                    std::string s_of_kE_lj_filename1 = 
                        output_dir + np_strings[nu] + "_s_of_k_" + 
                        "E_at_" + e_strm.str() + "MeV_" + L_array[L] 
                        + j_string + "2.out";

                    std::ofstream file1( s_of_kE_lj_filename1.c_str() );
*/
                    // Propagator
                    cmatrix_t G = B.propagator( rmesh, E, L, xj );

                    // Spectral Function in momentum space, S( k; E ) 
                    for( unsigned int nk = 0; nk < kmesh1.size(); ++nk ) {

                        double rsums = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                            double jl1 = bess_mtx1( nk, i );

                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            double jl2 = bess_mtx1( nk, j );
                    
                            rsums -= rmesh[i] * jl1 * imag( G( i, j ) ) 
                                   * rmesh[j] * jl2 * rdelt * 2 / M_PI / M_PI;
                        }
                        } // end loop over radial coordinates

//                        file1 << kmesh1[nk] << " " << rsums << std::endl;

                        S_of_kE_mtx1( nk, m ) += ( 2 * xj + 1 ) * rsums;

                    } // end loop over k

//                    file1.close();
//                    file1.clear();

                    // r-space spectral functions

                } // end loop over energy

                // Calculate S( k; E ) for E-slice 
                matrix_t S_of_kE_mtx_lj( kmesh2.size(), emesh2.size() );
                for ( unsigned int m = 0; m < emesh2.size(); ++m ) {
            
                    double E = emesh2[m];

                    // Propagator
                    cmatrix_t G = B.propagator( rmesh, E, L, xj );

                    // Spectral Function in momentum space, S( k; E ) 
                    for( unsigned int nk = 0; nk < kmesh2.size(); ++nk ) {

                        double rsums = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                            double jl1 = bess_mtx2( nk, i );

                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            double jl2 = bess_mtx2( nk, j );
                    
                            rsums -= rmesh[i] * jl1 * imag( G( i, j ) ) 
                                   * rmesh[j] * jl2 * rdelt * 2 / M_PI / M_PI;
                        }
                        } // end loop over radial coordinates

                        S_of_kE_mtx_lj( nk, m ) = ( 2 * xj + 1 ) * rsums;
                        S_of_kE_mtx2( nk, m ) += ( 2 * xj + 1 ) * rsums;

                    } // end loop over k

                } // end loop over energy

                S_of_kE_mtx2_lj_vec.push_back( S_of_kE_mtx_lj );

            }// end loop over j

        } // end loop over L

        // write out results summed over the lj combinations
        for ( unsigned int j = 0; j < emesh1.size(); ++j ) {

            std::ostringstream e_strm;
            e_strm << "m" << std::abs( emesh1[j] );
            std::string s_of_kE_filename1 = 
                output_dir + np_strings[nu] + "_s_of_k_" + 
                "E_at_" + e_strm.str() + "MeV.out"; 

            std::ofstream file3( s_of_kE_filename1.c_str() );

            for ( unsigned int i = 0; i < kmesh1.size(); ++i ) {

                file3 << kmesh1[i] << " " << S_of_kE_mtx1( i, j ) << std::endl;
            }

            file3.close();
            file3.clear();
        }

        // write in units of MeV^-4 sr^-1
        double fac = std::pow( hbarc, 3 ) * 4 * M_PI;
        for ( unsigned int i = 0; i < kmesh2.size(); ++i ) {
            std::ostringstream k_strm;
            k_strm << pmesh[i];
            std::string s_of_kE_filename2 = 
                output_dir + np_strings[nu] + "_s_of_E_" + 
                "p_at_" + k_strm.str() + "MeV_over_c.out"; 

            std::ofstream file4( s_of_kE_filename2.c_str() );

            for ( unsigned int j = 0; j < emesh2.size(); ++j ) {

                // write energy in terms of missing energy
                // and in units of GeV
                file4 << std::abs( emesh2[j] / 1000 ) << " " 
                      << S_of_kE_mtx2( i, j ) / fac << " "
                      << S_of_kE_mtx2( i, j ) / fac / calc_norm << std::endl;
            }

            file4.close();
            file4.clear();
        }

                // write in units of MeV^-4 sr^-1
        for ( unsigned int nk = 0; nk < kmesh2.size(); ++nk ) {

            for ( int L = 0; L < lmax + 1; ++L ) {
        
                    std::ostringstream k_strm;
                    k_strm << pmesh[nk];

                    std::string s_of_kE_lj_filename2 = 
                        output_dir + np_strings[nu] + "_s_of_E_p_at_" 
                        + k_strm.str() + "MeV_over_c_" + L_array[L] 
                        + ".out";

                    std::ofstream file2( s_of_kE_lj_filename2.c_str() );
                    for ( unsigned int m = 0; m < emesh2.size(); ++m ) {

                        // write energy in terms of missing energy
                        // and in units of GeV
                        if ( L == 0 ) {
                            file2 << std::abs( emesh2[m] / 1000 ) << " " 
                                  << S_of_kE_mtx2_lj_vec[0]( nk, m ) / fac 
                                  << std::endl;
                        }
                        else {

                            double jg = S_of_kE_mtx2_lj_vec[2*L]( nk, m );
                            double jl = S_of_kE_mtx2_lj_vec[2*L- 1]( nk, m );

                            file2 << std::abs( emesh2[m] / 1000 ) << " " 
                                  << jg / fac << " " << jl / fac << " "
                                  << ( jg + jl ) / fac << std::endl;
                        }
                    }

                    file2.close();
                    file2.clear();

            } // end loop over L

        } // end loop over momentum

    } // end loop over tz

    return 1;
}
