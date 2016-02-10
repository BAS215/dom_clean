#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include "read_parameters.h"
#include "legendre.h"
#include "boundRspace.h"
#include "meshes.h"
#include "density.h"
#include "io.h"
#include "Gauss_q.h"

int main() {

    double hbarc = 197.327; // hbar * c in [MeV fm]
    double Mass = 938; // nucleon mass in MeV

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

  //  std::string output_dir = "Output_" + parameters_string + "/Ener/";
  //  std::string parameters_filename = parameters_string + ".inp";

    std::string output_dir = "Output_test/S/";
   // std::string parameters_filename = "Output_hosfit/roott/Nov26/hosfitupdated.inp";
    std::string parameters_filename = "Input.inp";

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

    std::string L_array[8] = { "s", "p", "d", "f", "g", "h", "i", "j" };
    std::string J_array[8] = { "1", "3", "5", "7", "9", "11", "13", "15" };

    // Create momentum space grid
    std::vector<double> kmesh;
    std::vector<double> kweights;
    double const kmax = 6.0;
    int const kpts = 104;
    kmesh.resize( kpts );
    kweights.resize( kpts );
    GausLeg( 0., kmax, kmesh, kweights );

    // Create radial grid
    std::vector<double> rmesh;
    std::vector<double> rweights;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        rweights.push_back( rdelt );
    }

    // Prepare stuff for output files
    std::vector< std::string > np_strings;
    np_strings.push_back( n_string );
    np_strings.push_back( p_string );

    double total_energy = 0;
    double calc_A = 0;

    //Potential specifications
    int type = 1; // 1 is P.B. form (average), 0 is V.N. form
    int mvolume = 4;
    int AsyVolume = 1;

    // --- CALCULATIONS --- //

        std::string lj_particleN_share = output_dir + "ParticleN_shares.out";
        std::ofstream particleNumber_file( lj_particleN_share.c_str() );
    // Loop over protons and neutrons
    double Zp; // number of protons in projectile

    for ( unsigned int nu = 0; nu < Nu_vec.size(); ++nu ) {
        
        std::map< std::string, int > lj_map = map_vec[nu];

        double tz = nu - 0.5; // +0.5 for protons, -0.5 for neutrons
        if ( tz < 0 ) {Zp = 0;} else { Zp = 1;}
         particleNumber_file << tz << std::endl;

        // Nucleus Object
        const NuclearParameters &Nu = Nu_vec[nu];

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

        // Construct Potential Object
        pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
        pot *U1 = &U;
       // lmax = 0;
        // Construct boundRspace Object
//      double Emax = Nu.Ef - Nu.Wgap * p.fGap_B; // energy < Ef where imaginary part begins
        double Emax =  Nu.Ef;
        double Emin = -200 + Emax;
        boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, lmax , Nu.Z, Zp, Nu.A, U1 );

        double tol = 0.01;
        std::vector< lj_eigen_t > bound_levels = B.get_bound_levels( rmesh, tol );

        std::vector< mesh_t > emesh_vec = B.get_emeshes( rmesh, Emin, Emax, bound_levels );

        std::vector<double> n_of_k;
        n_of_k.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_qh;
        n_of_k_qh.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_c_tot;
        n_of_k_c_tot.assign( kmesh.size(), 0 );
        std::vector<double> n_of_k_qh_peak_tot;
        n_of_k_qh_peak_tot.assign( kmesh.size(), 0 );

        // Loop over lj channels
        for ( int l = 0; l < lmax + 1 ; ++l ) {

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

                std::vector<double> n_of_k_c_lj;
                std::vector<double> n_of_k_qh_lj;
                std::vector<double> n_of_k_qh_peak;
                std::vector<std::vector<double> > n_of_k_qh_peak_Matrice;
                n_of_k_c_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_lj.assign( kmesh.size(), 0 );
                n_of_k_qh_peak.assign( kmesh.size(), 0 );

	        matrix_t Slj(kmesh.size(),kmesh.size());
		for (int i=0;i<kmesh.size();++i){for (int j=0;j<kmesh.size();++j){Slj(i,j) = 0.;}}
	        matrix_t n_r_lj(rmesh.size(),rmesh.size());
		for (int i=0;i<rmesh.size();++i){for (int j=0;j<rmesh.size();++j){n_r_lj(i,j) = 0.;}}

                // Find self-consistent solutions
                 std::cout << "lj = " << l << " " << xj << std::endl;
                 std::vector<cmatrix_t> G;  
		 double E;
		 vector<double> Edelt;
                for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                    
                     E = emesh[m].first;
       		     Edelt.push_back(emesh[m].second);
                    // Propagator
	            G.push_back( B.propagator( rmesh, E, l, xj));
	        }
                // Loop over energy
		    for (int i = 0 ; i < rmesh.size(); ++i){
		         for (int j = 0 ; j < rmesh.size(); ++j){
			     double sumnr = 0.;
                             for ( unsigned int m = 0; m < emesh.size(); ++m ) {
			            sumnr += rmesh[i] * rmesh[j] * imag(G[m](i,j)) * Edelt[m];
		            }
		    n_r_lj(i,j) = sumnr;	
		         }
                    }
                  //  for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {
                    for( unsigned int mk = 0; mk < kmesh.size(); ++mk ) {
                         double rsums = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                            double jl1 = bess_mtx( mk, i );
                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            double jl2 = bess_mtx( mk, j );
                            
                            rsums -= n_r_lj(i,j) * jl1 * jl2 * rdelt;
                        }
                        } // end loop over radial coordinates
		//	Slj(nk,mk) = 2./M_PI/M_PI * rsums * ( 2 * xj + 1 )/( 4 * M_PI ); 		
			 n_of_k_c_lj[mk] += 2./M_PI/M_PI * rsums * ( 2 * xj + 1 )/( 4 * M_PI ); 		
		    }
                   // }    // integral over k
                     
                      // end loop over k
                 // end loop over energy
//                for (int nk = 0 ; nk < kmesh.size() ; ++nk){
//		     n_of_k_c_lj[nk] = Slj(nk,nk);} 
                // loop over the orbitals and write out the 
                // bound state information to a file
                for ( unsigned int N = 0; N < bound_info.size(); ++N ) {
                    // Quasiparticle energy
                    double QPE = bound_info[N].first; 

                    // Quasiparticle wavefunction
                    std::vector<double> &QPF = bound_info[N].second;

                    // Transform wavefunction to momentum space
                    std::vector<double> kQPF;
                    for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {

                        double sum = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {

                           sum += std::sqrt( 2 / M_PI ) * rweights[i] 
                                * std::pow( rmesh[i], 2 ) * QPF[i] * bess_mtx( nk, i );
                        }
                        kQPF.push_back( sum );
                      }

                    // This should automatically be normalized to N or Z
                    // since the wavefunctions are normalized to 1
                    if ( QPE < Nu.Ef ) {
                        for ( unsigned int kk = 0; kk < kmesh.size(); ++kk ) {

                            n_of_k_qh_lj[kk] += ( 2 * xj + 1 ) / ( 4 * M_PI ) * kQPF[kk] * kQPF[kk]; 
                        } // end loop over k
                    } 
                } // end loop over N 
                // Write Out Momentum Distribution for each lj channel
                // Also, add up the contribution from each channel to get 
                // the total momentum distribution
                std::string n_of_k_lj_filename = output_dir + "n_of_k_" 
                                               + np_strings[nu] 
                                               + "_" + L_array[l] + j_string 
                                               + "2.out";
                std::ofstream n_of_k_lj_file( n_of_k_lj_filename.c_str() );
		double NOFK_lj = 0;
                for ( unsigned int i = 0; i < kmesh.size(); ++i ) {
                    double n_of_k_lj = n_of_k_c_lj[i] + n_of_k_qh_peak[i];

                    n_of_k_lj_file << kmesh[i] << " " << n_of_k_lj 
                                   << " " << n_of_k_qh_lj[i] << std::endl;

                    n_of_k[i] += n_of_k_lj;
                    n_of_k_qh[i] += n_of_k_qh_lj[i];
                    n_of_k_qh_peak_tot[i] += n_of_k_qh_peak[i];
                    n_of_k_c_tot[i] += n_of_k_c_lj[i];

		    //Particle contribution to each state
		    // 
		    NOFK_lj += 4 * M_PI * kweights[i] * std::pow( kmesh[i], 2 ) * n_of_k_lj;
                }

                particleNumber_file <<index <<" " << "l= "<< l << " " <<"j= " << xj <<" " << NOFK_lj   << std::endl; 
				  //  <<" " <<n_of_k_qh_peak_tot[sarekar-1] <<" " <<n_of_k_c_tot[sarekar-1] <<std::endl;

                n_of_k_lj_file.close();
                n_of_k_lj_file.clear();
	    }//up
         }//lj
        // Output Files
        std::string strength_filename = output_dir + np_strings[nu] 
                                      + "_k_strength.out";
        std::ofstream strength_file( strength_filename.c_str() );

        std::string n_of_k_filename = output_dir + np_strings[nu]
                                    + "_momentum_dist.out";
        std::ofstream n_of_k_file( n_of_k_filename.c_str() ); 

        // Calculate Particle Number and Kinetic Energy
        double norm = 0;
        for ( unsigned int n = 0; n < kmesh.size(); ++n ) {
            norm += 4 * M_PI * kweights[n] * std::pow( kmesh[n], 2 ) * n_of_k[n];
        }
        // Energetics
        double N_or_Z; // actual number of neutrons or protons
        if( tz > 0 ) { 
            N_or_Z = Nu.Z;
        }
        else {
            N_or_Z = Nu.N();
        }

        calc_A += norm;
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

std::cout<< "HELLO tz LOOP" <<std::endl;
    }//tz
particleNumber_file.close();
return 1;
}

