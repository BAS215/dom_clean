//
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include "legendre.h"
#include "manybody.h"
#include "meshes.h"
#include "numerical.h"
#include "io.h"
#include "Gauss_q.h"


int main( ) {

    double hbarc = 197.327; // hbar * c in [MeV fm]
    double Mass = 938; // nucleon mass in MeV

    std::string input_dir = "Input/"; //"Input/";

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

    // Knowing the expected number of bound states is useful when
    // looking for particle states that are near the continuum

    // These maps hold the number of bound states for each lj
    std::map<std::string, int> n_lj_map; // neutrons
    std::map<std::string, int> p_lj_map; // protons
    std::vector< std::map< std::string, int > > map_vec;
    for ( int n = 0; n < num_lj; ++n ) {

        std::string key;
        int n_num_bound, p_num_bound;
        config_file >> key >> n_num_bound >> p_num_bound;

        std::pair< std::map< std::string, int >::iterator, bool > n_ret;
        std::pair< std::map< std::string, int >::iterator, bool > p_ret;
        n_ret = n_lj_map.insert ( std::make_pair( key, n_num_bound ) );
        p_ret = p_lj_map.insert ( std::make_pair( key, p_num_bound ) );

        if ( n_ret.second == false ) {
            
            std::cout << "element " << key << " already exists "; 
            std::cout << "with a value of " << n_ret.first->second << std::endl;
        }

        if ( p_ret.second == false ) {
            
            std::cout << "element " << key << " already exists "; 
            std::cout << "with a value of " << p_ret.first->second << std::endl;
        }
    }

    map_vec.push_back( n_lj_map );
    map_vec.push_back( p_lj_map );

    config_file.close();
    config_file.clear();

    std::string output_dir = "Output_" + parameters_string + "/tensory/";
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

    // Create momentum space grid
    std::vector<double> kmesh;
    std::vector<double> kweights;
    double const kmax = 6.0;
    int const kpts = 104;
    kmesh.resize( kpts );
    kweights.resize( kpts );
    GausLeg( 0., kmax, kmesh, kweights );

    // Prepare stuff for output files
    std::vector< std::string > np_strings;
    np_strings.push_back( n_string );
    np_strings.push_back( p_string );

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
        double Elim;
        if ( tz < 0 ) {

            Zp = 0;
            Elim = -62; // below 0s1/2 level for neutrons in 40Ca
        }
        else { 
            Zp = 1;
            Elim = -56; // below 0s1/2 level for protons in 40Ca
        }


        std::map< std::string, int > lj_map = map_vec[nu];

        // Nucleus Object
        const NuclearParameters &Nu = Nu_vec[nu];

        std::cout << "Ef = " << Nu.Ef << std::endl;
        std::cout << "A = " << Nu.A << std::endl;
        std::cout << "Z = " << Nu.Z << std::endl;

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

        // Construct Potential Object
        pot U = get_bobs_pot( type, mvolume, AsyVolume, tz, Nu, p );
        // Create Energy Grid for calculating density matrix etc.
      //  double Emax = Nu.Ef - Nu.gap * p.fGap;
        double Emax =   Nu.Ef;
        double Emin = -200 + Emax;
        double Elowmax = -60;
        double Emedmax = Emax - 21.; // 3 * Nu.gap;
        double Edeltlow = 1;// 2;
        double Edeltmed = .1;//1;
        double Edelthigh = 0.005;

        std::vector<std::pair< double, double > > emesh = 
            energy_mesh_fast( Emin, Elowmax, Emedmax, Emax, 
                              Edeltlow, Edeltmed, Edelthigh );

       //  Elim = Emax;

      //    std::vector<std::pair< double, double > > emesh = 
      //      energy_mesh_slow( Emin, Emax, 0.05 );



        std::vector<double> n_of_k;
        n_of_k.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_qh;
        n_of_k_qh.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_Elim;
        n_of_k_Elim.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_c_tot;
        n_of_k_c_tot.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_qh_peak_tot;
        n_of_k_qh_peak_tot.assign( kmesh.size(), 0 );

        // Loop over lj channels
        for ( int l = 0; l < lmax + 1; ++l ) {

            // Create Bessel Function matrix in k and r
            matrix_t bess_mtx( kmesh.size(), rmesh.size() );
            for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {
            for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

                double rho = kmesh[nk] * rmesh[nr];

                bess_mtx( nk, nr ) = gsl_sf_bessel_jl( l, rho );
            }
            }

            for( int up = -1; up < 2; up+=2 ) {

                double xj = l + up / 2.0;
                int j_index = ( up - 1 ) / 2;
                if ( xj < 0 ) continue;

                std::string j_string;
                if ( l == 0 ) j_string = J_array[ l ];
                else j_string = J_array[ l + j_index ];

                std::string lj_string = L_array[l] + j_string + "2";


                std::vector<double> n_of_k_c_lj;
                std::vector<double> n_of_k_c_Elim_lj;
                std::vector<double> n_of_k_qh_lj;
                std::vector<double> n_of_k_qh_peak;
                n_of_k_c_lj.assign( kmesh.size(), 0 );
                n_of_k_c_Elim_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_peak.assign( kmesh.size(), 0 );

                // Find self-consistent solutions
                int nb = lj_map[lj_string]; // # of bound states
                double Estart = Nu.Ef;

                std::vector<eigen_t> bound_info = 
                    slfcn_eig( rmesh, Estart, nb, l, xj, Nu, U );
                std::cout << "lj = " << l << " " << xj << std::endl;

                // Loop over energy
                double Esum = 0;
                double Esum2 = 0; // sum for energies up to E = Elim
                for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                    
                    double E = emesh[m].first;
                    double Edelt = emesh[m].second;

                    // Propagator
                    cmatrix_t G = propagator( rmesh, E, l, xj, Nu, U );

                    // Spectral Function in momentum space    
                    // Integrate E * S( k; E ) over k and E for 
                    // calculation of total energy
                    double ksum = 0;
                    for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {

                        double rsums = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                            double jl1 = bess_mtx( nk, i );

                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            double jl2 = bess_mtx( nk, j );
                            
                            rsums -= rmesh[i] * jl1 * imag( G( i, j ) ) *
                                     rmesh[j] * jl2 * rdelt * 2 / M_PI / M_PI;
                        }
                        } // end loop over radial coordinates

                        n_of_k_c_lj[nk] += Edelt * rsums * ( 2 * xj + 1 ) 
                                         / ( 4 * M_PI );

                        if ( E < Elim ) {

                            n_of_k_c_Elim_lj[nk] += 
                                Edelt * rsums * ( 2 * xj + 1 ) / ( 4 * M_PI );

                        }

                        // integral over k
                        ksum += kweights[nk] * std::pow( kmesh[nk], 2 )
                              * rsums;
                    } // end loop over k

                    // integral over energy
                    Esum += E * ksum * Edelt;

                    if ( E < Elim ) Esum2 += E * ksum * Edelt;

                } // end loop over energy



                // loop over the orbitals and write out the 
                // bound state information to a file
                for ( unsigned int N = 0; N < bound_info.size(); ++N ) {

                    // Quasiparticle energy
                    double QPE = bound_info[N].first; 

                    // Quasiparticle wavefunction
                    std::vector<double> &QPF = bound_info[N].second;

                    // Spectroscopic Factor
                    double S = sfactor( rmesh, QPE, l, xj, QPF, U );

                    // Transform wavefunction to momentum space
                    std::vector<double> kQPF;
                    for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {

                        double sum = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {

                           sum += std::sqrt( 2 / M_PI ) * rweights[i] 
                                * std::pow( rmesh[i], 2 ) * QPF[i] 
                                * bess_mtx( nk, i );
                        }
                        kQPF.push_back( sum );
                    }

                    // This should automatically be normalized to N or Z
                    // since the wavefunctions are normalized to 1
                    if ( QPE < Nu.Ef ) {
                        for ( unsigned int kk = 0; kk < kmesh.size(); ++kk ) {

                            n_of_k_qh_lj[kk] += ( 2 * xj + 1 ) / ( 4 * M_PI )
                                              * kQPF[kk] * kQPF[kk]; 
                        } // end loop over k
                    }
        
                    // Peak contributions
                    if ( ( QPE > Emax ) && ( QPE < Nu.Ef ) ) {

                    //    T_plus_2V_qh += S * QPE;

                        // This needs to be added to the continuum part
                        for ( unsigned int kk = 0; kk < kmesh.size(); ++kk ) {

                            n_of_k_qh_peak[kk] += ( 2 * xj + 1 ) / ( 4 * M_PI )
                                                * S * kQPF[kk] * kQPF[kk];

                        } // end loop over k

                        // Add peak contributions to the momentum distribution
                        std::cout << "added to momentum distribution " 
                                  << N << " " << l << " " << xj << " " 
                                  << QPE << " " << S << std::endl;

                    } // endif

                } // end loop over N 




                // Write Out Momentum Distribution for each lj channel
                // Also, add up the contribution from each channel to get 
                // the total momentum distribution
                std::string n_of_k_lj_filename = output_dir + "n_of_k_" 
                                               + np_strings[nu] 
                                               + "_" + L_array[l] + j_string 
                                               + "2.out";
                std::ofstream n_of_k_lj_file( n_of_k_lj_filename.c_str() );

                for ( unsigned int i = 0; i < kmesh.size(); ++i ) {
                    double n_of_k_lj = n_of_k_c_lj[i] + n_of_k_qh_peak[i];

                    n_of_k_lj_file << kmesh[i] << " " << n_of_k_lj 
                                   << " " << n_of_k_qh_lj[i] << std::endl;

                    n_of_k[i] += n_of_k_lj;
                    n_of_k_qh[i] += n_of_k_qh_lj[i];
                    n_of_k_Elim[i] += n_of_k_c_Elim_lj[i];
                    n_of_k_qh_peak_tot[i] += n_of_k_qh_peak[i];
                    n_of_k_c_tot[i] += n_of_k_c_lj[i];


                }

                n_of_k_lj_file.close();
                n_of_k_lj_file.clear();


            } //end loop over j
        } // end loop over l

        // Output Files
        std::string strength_filename = output_dir + np_strings[nu] 
                                      + "_k_strength.out";
        std::ofstream strength_file( strength_filename.c_str() );

        std::string n_of_k_filename = output_dir + np_strings[nu]
                                    + "_momentum_dist.out";
        std::ofstream n_of_k_file( n_of_k_filename.c_str() ); 

        double N_or_Z; // actual number of neutrons or protons
        if( tz > 0 ) { 
            N_or_Z = Nu.Z;
        }
        else {
            N_or_Z = Nu.N();
        }

        double norm = 0;
        for ( unsigned int n = 0; n < kmesh.size(); ++n ) {
            norm += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 2 )
                  * n_of_k[n];
              } 



        // Momentum Distribution and strength
        double k_strength = 0;
        double qh_strength = 0;
        for ( unsigned int n = 0; n < kmesh.size(); ++n ) {
            n_of_k_file << kmesh[n] << " " << n_of_k[n] / norm        
                                    << " " << n_of_k_qh[n] / N_or_Z 
                                    << std::endl;

            //verify that n_of_k_qh is normalized to N or Z
            qh_strength += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 2 ) 
                     * n_of_k_qh[n] / N_or_Z;

            k_strength += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 2 ) 
                        * n_of_k[n] / norm;
            
            strength_file << kmesh[n] << " " << k_strength 
                                      << " " << qh_strength << std::endl;
        }

        strength_file << " " << std::endl;
        strength_file << "norm = " << norm << std::endl;
        strength_file << "N or Z = " << N_or_Z << std::endl;



        n_of_k_file.close();
        n_of_k_file.clear();
        strength_file.close();
        strength_file.clear();


    } // end loop over tz

    return 1;
}
