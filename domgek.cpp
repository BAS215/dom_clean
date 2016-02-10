//
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
//#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include "legendre.h"
#include "boundRspace.h"
//#include "manybody.h"
#include "boundRspace.h"
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

    std::string output_dir = "Output_" + parameters_string + "/Energy/";
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

    std::string energy_filename = output_dir + nucleus_string + "_energy.out";
    std::ofstream energy_file( energy_filename.c_str() );
    double total_energy = 0;
    double calc_A = 0;

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
        pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
        pot *U1 = &U;
        // Create Energy Grid for calculating density matrix etc.
        // Construct boundRspace Object
        double Emax = Nu.Ef - Nu.Wgap * p.fGap_B; // energy < Ef where imaginary part begins
        double Emin = -200 + Emax;
        boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, lmax, Nu.Z, Zp, Nu.A, U1 );

        double tol = 0.01;
        std::vector< lj_eigen_t > bound_levels = B.get_bound_levels( rmesh, tol );

        std::vector< mesh_t > emesh_vec = B.get_emeshes( rmesh, Emin, Emax, bound_levels );


        double T_plus_2V = 0; // for total energy calculation
        double T_plus_2V_qh_peak_tot = 0;
        double T_plus_2V_c_tot = 0;
        std::vector< double > T_plus_2V_c_lj_vec;
        std::vector< double > T_plus_2V_c_Elim_lj_vec;
        std::vector< double > e_kin_c_lj_vec;
        std::vector< double > e_kin_peak_lj_vec;
        std::vector< double > e_kin_c_Elim_lj_vec;

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

                int index = B.index_from_LJ( l, xj );
                mesh_t &emesh = emesh_vec[index];
                std::vector< eigen_t > &bound_info = bound_levels[index];

                std::string lj_string = L_array[l] + j_string + "2";

                double T_plus_2V_c = 0; // contribution from continuum
                double T_plus_2V_qh = 0; // contribution from quasiholes
                double T_plus_2V_c_Elim = 0; // contribution for E < Elim

                std::vector<double> n_of_k_c_lj;
                std::vector<double> n_of_k_c_Elim_lj;
                std::vector<double> n_of_k_qh_lj;
                std::vector<double> n_of_k_qh_peak;
                n_of_k_c_lj.assign( kmesh.size(), 0 );
                n_of_k_c_Elim_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_peak.assign( kmesh.size(), 0 );

                // Find self-consistent solutions
    //            int nb = lj_map[lj_string]; // # of bound states
                double Estart = Nu.Ef;

                std::cout << "lj = " << l << " " << xj << std::endl;

                // Loop over energy
                double Esum = 0;
                double Esum2 = 0; // sum for energies up to E = Elim
                for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                    
                    double E = emesh[m].first;
                    double Edelt = emesh[m].second;

                    // Propagator
                    cmatrix_t G = B.propagator( rmesh, E, l, xj );

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

                T_plus_2V_c = Esum;
                T_plus_2V_c_Elim = Esum2;

                T_plus_2V_c_lj_vec.push_back( ( 2 * xj + 1 ) * T_plus_2V_c );
                T_plus_2V_c_Elim_lj_vec.push_back( 
                   ( 2 * xj + 1 ) * T_plus_2V_c_Elim );


                // loop over the orbitals and write out the 
                // bound state information to a file
                for ( unsigned int N = 0; N < bound_info.size(); ++N ) {

                    // Quasiparticle energy
                    double QPE = bound_info[N].first; 

                    // Quasiparticle wavefunction
                    std::vector<double> &QPF = bound_info[N].second;

                    // Spectroscopic Factor
                    double S = B.sfactor( rmesh, QPE, l, xj, QPF );

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

                       T_plus_2V_qh += S * QPE;

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

                T_plus_2V += ( 2 * xj + 1 ) * ( T_plus_2V_c + T_plus_2V_qh );

                T_plus_2V_qh_peak_tot += ( 2 * xj + 1 ) * T_plus_2V_qh;
                T_plus_2V_c_tot += ( 2 * xj + 1 ) * T_plus_2V_c;


                // Write Out Momentum Distribution for each lj channel
                // Also, add up the contribution from each channel to get 
                // the total momentum distribution
                std::string n_of_k_lj_filename = output_dir + "n_of_k_" 
                                               + np_strings[nu] 
                                               + "_" + L_array[l] + j_string 
                                               + "2.out";
                std::ofstream n_of_k_lj_file( n_of_k_lj_filename.c_str() );

                double e_kin_c_lj = 0;
                double e_kin_c_Elim_lj = 0;
                double e_kin_peak_lj = 0;
                for ( unsigned int i = 0; i < kmesh.size(); ++i ) {
                    double n_of_k_lj = n_of_k_c_lj[i] + n_of_k_qh_peak[i];

                    n_of_k_lj_file << kmesh[i] << " " << n_of_k_lj 
                                   << " " << n_of_k_qh_lj[i] << std::endl;

                    n_of_k[i] += n_of_k_lj;
                    n_of_k_qh[i] += n_of_k_qh_lj[i];
                    n_of_k_Elim[i] += n_of_k_c_Elim_lj[i];
                    n_of_k_qh_peak_tot[i] += n_of_k_qh_peak[i];
                    n_of_k_c_tot[i] += n_of_k_c_lj[i];

                    e_kin_c_lj += kweights[i] * std::pow( kmesh[i], 4 ) 
                                * std::pow( hbarc, 2 ) / ( 2 * Mass ) 
                                * n_of_k_c_lj[i] * 4 * M_PI; 

                    e_kin_peak_lj += kweights[i] * std::pow( kmesh[i], 4 ) 
                                   * std::pow( hbarc, 2 ) / ( 2 * Mass ) 
                                   * n_of_k_qh_peak[i] * 4 * M_PI; 

                    e_kin_c_Elim_lj += kweights[i] * std::pow( kmesh[i], 4 ) 
                                     * std::pow( hbarc, 2 ) / ( 2 * Mass ) 
                                     * n_of_k_c_Elim_lj[i] * 4 * M_PI;

                }

                n_of_k_lj_file.close();
                n_of_k_lj_file.clear();

                e_kin_c_lj_vec.push_back( e_kin_c_lj );
                e_kin_peak_lj_vec.push_back( e_kin_peak_lj );
                e_kin_c_Elim_lj_vec.push_back( e_kin_c_Elim_lj );

            } //end loop over j
        } // end loop over l

        // Output Files
        std::string strength_filename = output_dir + np_strings[nu] 
                                      + "_k_strength.out";
        std::ofstream strength_file( strength_filename.c_str() );

        std::string n_of_k_filename = output_dir + np_strings[nu]
                                    + "_momentum_dist.out";
        std::ofstream n_of_k_file( n_of_k_filename.c_str() ); 

        // Calculate Particle Number and Kinetic Energy
        double norm = 0;
        double kineticE = 0;
        double kineticE_c = 0;
        double kineticE_peak = 0;
        for ( unsigned int n = 0; n < kmesh.size(); ++n ) {
            norm += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 2 ) 
                  * n_of_k[n];

            kineticE += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 4 ) 
                      * n_of_k[n] * std::pow( hbarc, 2 ) / ( 2 * Mass );

            kineticE_c += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 4 ) 
                        * n_of_k_c_tot[n] * std::pow( hbarc, 2 ) / ( 2 * Mass );

            kineticE_peak += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 4 ) 
                           * n_of_k_qh_peak_tot[n] 
                           * std::pow( hbarc, 2 ) / ( 2 * Mass );

        }

        // Energetics
        double E_total = ( kineticE + T_plus_2V ) / 2;
        double E_potential = ( T_plus_2V - kineticE ) / 2;

        double E_total_c4 = ( kineticE_c + T_plus_2V_c_tot ) / 2;
        double E_total_peak = ( kineticE_peak + T_plus_2V_qh_peak_tot ) / 2;

        if ( tz > 0 ) energy_file << "Protons " << std::endl;
        else energy_file << "Neutrons " << std::endl;

        energy_file << " " << std::endl;
        energy_file << "Total Energy: " << E_total << std::endl;
        energy_file << "Kinetic Energy: " << kineticE << std::endl;
        energy_file << "Potential Energy: " << E_potential << std::endl;
        energy_file << "T_plus_2V: " << T_plus_2V << std::endl;
        energy_file << " " << std::endl;
        energy_file << "Kinetic Energy from continuum: " << kineticE_c 
                    << std::endl;
        energy_file << "Kinetic Energy from qh delta peaks: " << kineticE_peak 
                    << std::endl;
        energy_file << "T_plus_2V from continuum: " << T_plus_2V_c_tot 
                    << std::endl;
        energy_file << "T_plus_2V from qh delta peaks: " << T_plus_2V_qh_peak_tot 
                    << std::endl;
        energy_file << "Total energy from continuum: " << E_total_c4 
                    << std::endl;
        energy_file << "Total energy from qh delta peaks: " << E_total_peak 
                    << std::endl;

        energy_file << " " << std::endl;
        double N_or_Z; // actual number of neutrons or protons
        if( tz > 0 ) { 
            energy_file << "E / Z = " << E_total / norm << std::endl;
            energy_file << "Z = " << norm << std::endl;
            N_or_Z = Nu.Z;
        }
        else {
            energy_file << "E / N = " << E_total / norm << std::endl;
            energy_file << "N = " << norm << std::endl;
            N_or_Z = Nu.N();
        }

        total_energy += E_total;
        calc_A += norm;

        energy_file << " " << std::endl;
        energy_file << "// --- Continuum Contributions to energy --- //" 
                    << std::endl;
        energy_file << "L  j  KE  KE [-inf,-50]  Total E  Total E [-inf,-50]"
                    << std::endl;

        
        double E_total_c = 0;
        double E_total_c2 = 0;
        double E_total_c3 = 0;
        double kineticE2_check = 0;

        for ( int L = 0; L < lmax + 1; ++L ) {
        for ( int nj = -1; nj <= 0; ++nj ) {
            
            double xj = L + 0.5 + nj;
            int lj_index = 2 * L + nj;
            if ( lj_index < 0 ) continue;

            double ekin = e_kin_c_lj_vec.at( lj_index );
            double ekin_Elim = e_kin_c_Elim_lj_vec.at( lj_index );
            double T_plus_2V_c = T_plus_2V_c_lj_vec.at( lj_index );
            double T_plus_2V_c_Elim = T_plus_2V_c_Elim_lj_vec.at( lj_index );

            double total = ( ekin + T_plus_2V_c ) / 2.0;
            double total_Elim = ( ekin_Elim + T_plus_2V_c_Elim ) / 2.0;
            energy_file << L << " " << xj << " " << ekin << " " << ekin_Elim 
                        << " " << total << " " << total_Elim << std::endl;

            E_total_c += total; // Should be the same as Etotal
            kineticE2_check += ekin_Elim;

            if ( L == 0 ) E_total_c2 += total_Elim;
            else E_total_c2 += total;

            E_total_c3 += total_Elim;
        }
        }

        energy_file << "Total Continuum Contribution = " 
                    << E_total_c << ", " << E_total_c / E_total * 100 
                    << "\%" << std::endl;

        energy_file << "Total Continuum Contribution ( s wave, " 
                    << "Emax = " << Elim << " ) = " << E_total_c2 
                    << ", " << E_total_c2 / E_total * 100 
                    << "\%" << std::endl;

        energy_file << "Total Continuum Contribution ( " 
                    << "Emax = " << Elim << " ) = " << E_total_c3 
                    << ", " << E_total_c3 / E_total * 100 
                    << "\%" << std::endl;

        energy_file << " " << std::endl;
        energy_file << "-----------------------------------------" << std::endl;
        energy_file << " " << std::endl;

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

    energy_file << " " << std::endl;
    energy_file << "Total Binding Energy: " << - total_energy << std::endl;
    energy_file << "E / A = " << total_energy / calc_A << std::endl;

    energy_file.close();

    return 1;
}
