//
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


int main() {

    double wightt_J [11] = { 2., 4., 2., 4., 6., 6., 8., 8., 10., 10., 12. };
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

    int num_so;
    config_file >> num_so;

    // These maps determine which orbits are included in the
    // spin-orbit correction to the charge density
    std::map<std::string, int> n_so_map; // neutrons
    std::map<std::string, int> p_so_map; // protons
    std::vector< std::map< std::string, int > > so_map_vec;
    for ( int n = 0; n < num_so; ++n ) {

        std::string key;

        // 1 if shell is filled, 0 if not at all (in IPM),
        // in between if shell is partially filled
        int n_shell, p_shell; 
        config_file >> key >> n_shell >> p_shell;

        std::pair< std::map< std::string, int >::iterator, bool > n_ret;
        std::pair< std::map< std::string, int >::iterator, bool > p_ret;
        n_ret = n_so_map.insert ( std::make_pair( key, n_shell ) );
        p_ret = p_so_map.insert ( std::make_pair( key, p_shell ) );

        if ( n_ret.second == false ) {
            
            std::cout << "element " << key << " already exists "; 
            std::cout << "with a value of " << n_ret.first->second << std::endl;
        }

        if ( p_ret.second == false ) {
            
            std::cout << "element " << key << " already exists "; 
            std::cout << "with a value of " << p_ret.first->second << std::endl;
        }
    }

    so_map_vec.push_back( n_so_map );
    so_map_vec.push_back( p_so_map );

    config_file.close();
    config_file.clear();

    bool hole;
    if ( fit_ph == 0 ) hole = true;
    else hole = false;

    std::string output_dir = "Output_" + parameters_string + "/Output_domgen/";
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

    std::string L_array[8] = { "s", "p", "d", "f", "g", "h", "i", "j" };
    std::string J_array[8] = { "1", "3", "5", "7", "9", "11", "13", "15" };


    // Create radial grid
    std::vector<double> rmesh;
    std::vector<double> rweights;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        rweights.push_back( rdelt );
    }

    // initialize charge density vector
    std::vector<double> chd_ls; //relativistic spin-orbit correction
    std::vector<double> chd_ls_test; //relativistic spin-orbit correction
    chd_ls.assign( rmesh.size(), 0 ); 
    chd_ls_test.assign( rmesh.size(), 0 ); 
   

    double E_c = -20;
    double HOSS[20][rmesh.size()];
//    double HOSS1[20][rmesh.size()];    
    for (unsigned int ii = 0 ;ii<20 ; ++ii){for (unsigned int k = 0 ; k<rmesh.size(); ++k) {HOSS[ii][k]=0.0;}}
    std::vector<vector < double > > HOSS1;   
    std::vector< double > EHOSS1; 
    

    std::string chd_filename = 
        output_dir + nucleus_string + "_charge_density.out";
    std::ofstream chd_file( chd_filename.c_str() );
    std::vector< std::vector<double> > chdf_np_vec;

    // Prepare stuff for output files
    std::vector< std::string > np_strings;
    np_strings.push_back( n_string );
    np_strings.push_back( p_string );

    //Potential specifications
    int type = 1; // 1 is P.B. form (average), 0 is V.N. form
    int mvolume = 4;
    int AsyVolume = 1;

    // --- CALCULATIONS --- //

    // Loop over protons and neutrons
    double Zp; // number of protons in projectile
    std::vector< double > coulomb_vec;
    coulomb_vec.assign( rmesh.size(), 0 ); // store coulomb potential
    for ( unsigned int nu = 0; nu < Nu_vec.size(); ++nu ) {
        
        std::map< std::string, int > lj_map = map_vec[nu];
        std::map< std::string, int > so_map = so_map_vec[nu];

        double tz = nu - 0.5; // +0.5 for protons, -0.5 for neutrons
        double mag_moment;

        if ( tz < 0 ) {
            Zp = 0;
            mag_moment = -1.91; // mu_n;
        }
        else {
            Zp = 1;
            mag_moment = 2.29; // Actually mu_p - 0.5

        }

        // Nucleus Object
        const NuclearParameters &Nu = Nu_vec[nu];

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

        // Construct Potential Object
        pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
        pot *U1 = &U;

        // Construct boundRspace Object
        double Emax = Nu.Ef - Nu.Wgap * p.fGap_B; // energy < Ef where imaginary part begins
        double Emin = -200 + Emax;
        boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, lmax, Nu.Z, Zp, Nu.A, U1 );

        double tol = 0.01;
        std::vector< lj_eigen_t > bound_levels = 
            B.get_bound_levels( rmesh, tol );

        std::vector< mesh_t > emesh_vec = 
            B.get_emeshes( rmesh, Emin, Emax, bound_levels );

        // store coulomb potential in order to compare with
        // proton charge distribution
        if ( tz > 0 ) {
            for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                coulomb_vec[i] += U.coulomb( rmesh[i] );
            }
        }

        // Vector for holding neutron or proton point distribution
        std::vector<double> point_dist;
        std::vector<double> point_dist_wfx; // calculated using wavefunctions
        point_dist.assign( rmesh.size(), 0 );
        point_dist_wfx.assign( rmesh.size(), 0 );

        // Output Files
        std::string qp_filename = 
            output_dir + np_strings[nu] + "_quasiparticle" + ".out";
        std::ofstream qp_file( qp_filename.c_str() );

        std::string strength_filename = 
            output_dir + np_strings[nu] + "_strengths.out";
        std::ofstream strength_file( strength_filename.c_str() );

        std::string density_filename = 
            output_dir + np_strings[nu] + "_density_dist.out";
        std::ofstream density_file( density_filename.c_str() ); 

        


        // Loop over lj channels
        double rsum = 0; // running sum of strength
        for ( int l = 0; l < lmax + 1; ++l ) {

            for( int up = -1; up < 2; up+=2 ) {

                double xj = l + up / 2.0;
                int j_index = ( up - 1 ) / 2;
                if ( xj < 0 ) continue;

                int index = B.index_from_LJ( l, xj );
                mesh_t &emesh = emesh_vec[index];
                std::vector< eigen_t > &bound_info = bound_levels[index];

                double lp; 
                if ( xj > l ) lp = l;
                else lp = - l - 1;
                
                std::string j_string;
                if ( l == 0 ) j_string = J_array[ l ];
                else j_string = J_array[ l + j_index ];

                std::string lj_string = L_array[l] + j_string + "2";

                std::string spf_filename;
                std::string chd_lj_filename;

                spf_filename = output_dir + np_strings[nu] + "_spectral_f_" 
                             + L_array[l] + j_string + "2.out";

                chd_lj_filename = output_dir + np_strings[nu] + "_chd_" 
                                + L_array[l] + j_string + "2.out";

                std::ofstream spf_file( spf_filename.c_str() );
                std::ofstream chd_lj_file( chd_lj_filename.c_str() );

                // Find self-consistent solutions
                std::cout << "lj = " << l << " " << xj << std::endl;




		std::vector<double> sumG;
                sumG.assign( rmesh.size() ,0);
               
               



                 
                // Loop over energy
                double strength = 0;
                matrix_t d_mtx( rmesh.size(), rmesh.size() ); // density matrix
                d_mtx.clear();
                std::vector<cmatrix_t> G_vec;
                for ( unsigned int m = 0; m < emesh.size(); ++m ) {
                    
                    double E = emesh[m].first;
                    double Edelt = emesh[m].second;

                    // Propagator
                    cmatrix_t G = B.propagator( rmesh, E, l, xj );

                    G_vec.push_back( G );

                    // Spectral Function 
                    complex_t trace = 0;
                    for( unsigned int i = 0; i < rmesh.size(); ++i ) {
				
				if (tz>0){ 
                                	sumG[i] = sumG[i] - 1/M_PI *imag( G(i,i))/(rmesh[i] * rmesh [i] * rdelt) * Edelt;
				}

                       
                        trace += G( i, i ); 
                    }

		    if( (tz>0) && (index==4) && (E>E_c)) {
							  HOSS1.push_back(sumG);
							  EHOSS1.push_back(E);
							  E_c+=.2;}			

                    double spf = -imag( trace ) / M_PI; 
                    spf_file << E << " " << spf << std::endl;

                    // Density Matrix
                    for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {
                            
                            //G( i, j ) /= rmesh[i] * rmesh[j] * rdelt;
                            d_mtx( i, j ) -= Edelt * imag( G( i, j ) ) / M_PI
                                           / ( rmesh[i] * rmesh[j] * rdelt );
                        }
                    }

                    strength += Edelt * spf;

                } // end loop over energy

          	if (tz > 0){
                    for (unsigned int i = 0 ; i<rmesh.size() ; ++i){
                        HOSS[index][i] = sumG[i];   
		
                    }
                } 
                // Charge Density
                std::vector<double> chd_lj;
                for( unsigned int i = 0; i < rmesh.size(); ++i ) {
    
                    double chd = ( 2 * xj + 1 ) * d_mtx( i, i ) / ( 4 * M_PI );
                    chd_lj.push_back( chd );
                }

                // loop over the orbitals and write out the 
                // bound state information to a file
                int intj = static_cast<int>( 2 * xj );
                for ( int N = 0; N < bound_info.size(); ++N ) {

                    std::string level_str = util::IntToStr( N ) + L_array[l] 
                                          + j_string + "2";
/*
                    std::string wfx_filename = output_dir + np_strings[nu] 
                        + "_" + level_str + ".wfx";
                    std::ofstream wfx_file( wfx_filename.c_str() );
*/
                    std::string chd_nlj_filename = output_dir + np_strings[nu] 
                        + "_" + level_str + ".chd";
                    std::ofstream chd_nlj_file ( chd_nlj_filename.c_str() );

                    std::string spf_nlj_filename = output_dir + np_strings[nu] 
                        + "_" + level_str + ".spf"; 
                    std::ofstream spf_nlj_file( spf_nlj_filename.c_str() );
                
                    // Quasiparticle energy
                    double QPE = bound_info[N].first; 

                    // Quasiparticle wavefunction
                    std::vector<double> &QPF = bound_info[N].second;

                    // Spectroscopic Factor
                    double S = B.sfactor( rmesh, QPE, l, xj, QPF );
                    double radius = B.rms_radius( rmesh, QPF );
        
                    // Occupation of continuum part
                    double occ = B.occupation( rmesh, d_mtx, QPF );

                    // fold spectral function with quasihole wavefunctions
                    double occ2 = 0; // should be the same as occ
                    //matrix_t projected_dmtx( rmesh.size(), rmesh.size() );
                    for( unsigned int g = 0; g < G_vec.size(); ++g ) {

                        double E = emesh[g].first;
                        double Edelt = emesh[g].second;

                        double fspf = 0;
                        for( unsigned int i = 0; i < rmesh.size(); ++i ) {
                        for( unsigned int j = 0; j < rmesh.size(); ++j ) {

                           fspf -= rmesh[i] * rmesh[j] * rdelt / M_PI
                                 * QPF[i] * QPF[j] * imag( G_vec[g]( i, j ) );
                        }
                        }

                        occ2 += fspf * Edelt;
                        spf_nlj_file << E << " " << fspf << std::endl;
                    }

                    // get rough estimate of quasihole width
                    double gamma = B.approx_width( rmesh, QPE, l, xj, QPF );

                    spf_nlj_file.close();
                    spf_nlj_file.clear();

                    qp_file << N << L_array[l] << intj << "/2" 
                            << " " << QPE << " " << S 
                            << " " << radius << " " << occ 
                            << " " << gamma << std::endl;


                    /* CHARGE DENSITY STUFF */

                    // calculate charge density using wavefunctions
                    std::vector<double> chd_nlj;
                    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
                        
                        chd_nlj.push_back( ( 2 * xj + 1 ) / ( 4 * M_PI )
                                            * QPF[i] * QPF[i] );
                                            
                        point_dist_wfx[i] += chd_nlj[i];
                        // Write out wavefunctions to a file
                        //wfx_file << rmesh[i] << " " << QPF[i] << std::endl;
                        chd_nlj_file << rmesh[i] << " " << chd_nlj[i] << " " 
                                     << rmesh[i] * chd_nlj[i] << std::endl;
                    }

                    //wfx_file.close();
                    //wfx_file.clear();
                    chd_nlj_file.close();
                    chd_nlj_file.clear();

                    double total_occ;
                    if ( ( QPE > Emax ) && ( QPE <= ( Nu.Ef ) ) ) {

                        strength += S;
                        total_occ = occ + S;

                        for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                            chd_lj[i] += ( 2 * xj + 1 ) * S * QPF[i] * QPF[i] 
                                       / ( 4 * M_PI );

                        } // end loop over r

                        // Add peak contributions to the charge density
                        std::cout << "added to charge density " << N << " " 
                                  << l << " " << xj << " " << QPE << " " 
                                  << S << std::endl;

                    } // endif
                    else total_occ = occ;

                    // calculate derivative of r * charge density for 
                    // relativistic spin-orbit correction. Taken from
                    // J. Phys. G: Nucl. Phys. 5, 1655 ( 1979 )

                    if( ( so_map.count(level_str) > 0 ) && 
                        ( so_map[level_str] > 0 ) ) {
                        
                        std::cout << "doing spinorbit correction" << std::endl;
                    for ( unsigned int i = 1; i < rmesh.size()-1; ++i ) {

                        double der = 
                            der_3pt( rdelt, rmesh[i-1] * chd_nlj[i-1], 
                                     rmesh[i+1] * chd_nlj[i+1] );

                        // correction
                        chd_ls[i] += 0.5 * total_occ * lp * mag_moment
                                   * der / std::pow( rmesh[i], 2 )
                                   * so_map[level_str];
/*
                        double der4 = 
                            der_4pt( rdelt, rmesh[i-2] * chd_nlj[i-2],
                                     rmesh[i-1] * chd_nlj[i-1],
                                     rmesh[i+1] * chd_nlj[i+1],
                                     rmesh[i+2] * chd_nlj[i+2] );
                        chd_ls_test[i] += 0.5 * total_occ * lp * mag_moment
                                        * der4 / std::pow( rmesh[i], 2 );
                        
*/
                    }
                    } // end if so_correction

                } // end loop over N 

                rsum += strength * ( 2 * xj + 1 );

                strength_file << L_array[l] << intj << "/2"
                              << " " << strength << " " << rsum << std::endl;

                // print out point distribution for each lj
                for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                    chd_lj_file << rmesh[i] << " " << chd_lj[i] << std::endl;
                    point_dist[i] += chd_lj[i];
                }

                spf_file.close();
                spf_file.clear();
                chd_lj_file.close();
                chd_lj_file.clear();

            } //end loop over j
        } // end loop over l

        qp_file.close();
        qp_file.clear();

        // Print to file the matter distribution. Find Normalization first
        double point_norm = 0;
        for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

            point_norm += 4 * M_PI * point_dist[i] * std::pow( rmesh[i], 2 )
                        * rdelt;
        }
std::cout<<point_dist[0]<<std::endl;
        double particle_number;
        if ( tz > 0 ) particle_number = Nu.Z;
        else particle_number = Nu.N();

        double norm_fac = particle_number / point_norm;
        double new_norm = 0; 
        for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
            point_dist[i] *= norm_fac;
            new_norm += 4 * M_PI * std::pow( rmesh[i], 2 ) * rdelt 
                      * point_dist[i];
        }

        strength_file << " " << std::endl;
        strength_file << "Number of Particles from point distribution: " 
                      << point_norm << std::endl;

        strength_file << "Number after normalization: " << new_norm 
                      << std::endl;

        strength_file.close();
        strength_file.clear();


        std::cout << "Folding density distribution" << std::endl;
        // Fold point charge density with neutron and proton charge distribution

        std::vector<double> chdf = 
            folded_ch_density( rmesh, rweights, point_dist, tz, Nu.A );

        chdf_np_vec.push_back( chdf );

        // approximate finite size of neutron with proton charge
        // distribution in order to calculate the neutron matter
        // distribution
        std::vector<double> matter_dist =
            matter_distribution( rmesh, rweights, point_dist, Nu.A );

        // write out point distribution and matter distribution to file
        for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
            density_file << rmesh[i] << " " << point_dist[i] << " " 
                         << point_dist_wfx[i] << " " 
                         << matter_dist[i] << std::endl;
        }

        density_file.close();
        density_file.clear();

        if (tz>0) {
	std::string outputdirrr="/scratch/mmahzoon/Jafer/s";
        std::ofstream spectral_r_file( "/scratch/mmahzoon/Jafer/spec_file_r.out" );
        std::ofstream Espectral_s_r_file( "/scratch/mmahzoon/Jafer/Espec_file_s_r.out" );
        for (int i=0 ; i<rmesh.size();++i) {  		
							double tot_HOSS = 0;	
							for(int jj=0 ; jj <11 ; ++jj)  { tot_HOSS+= wightt_J[jj] * HOSS[jj][i]; }
							spectral_r_file <<rmesh[i] << " " << HOSS[0][i] << " " << HOSS[1][i] << " " << HOSS[2][i] 
								   << " " << HOSS[3][i] << " " << HOSS[4][i] << " " << HOSS[5][i]
								   << " " << HOSS[6][i] << " " << HOSS[7][i] << " " << HOSS[8][i] 
								   << " " << HOSS[9][i] << " " << HOSS[10][i] <<" " <<tot_HOSS<<std::endl;}



         for (int i=0 ; i<HOSS1.size();++i) { 
                Espectral_s_r_file << i <<" "<< EHOSS1[i] << std::endl;
		std::string Stringg = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
		std::string spectral_s_r_file= outputdirrr + Stringg + ".out";
		std::ofstream spectral_sr_file( spectral_s_r_file.c_str() );
   	 	for (int j=0 ; j<rmesh.size() ; ++j) { 
    			spectral_sr_file << rmesh[j]<<" " << HOSS1[i][j] <<  std::endl;
 
	   	}   
		spectral_sr_file.close();	

	}   
        Espectral_s_r_file.close();
	spectral_r_file.close();	}
    } // end loop over Nu_vec

    if ( chdf_np_vec.size() < 2 ) { 
        std::cout << "chdf_np_vec has wrong size." << std::endl;
        throw std::exception();
    }

    std::vector<double> chd_tot;
    chd_tot.assign( rmesh.size(), 0 );

    double hbar_over_mc = 0.21;
    for( unsigned int i = 0; i < rmesh.size(); ++i ) {

        // add folded proton charge distribution
        chd_tot[i] += chdf_np_vec[1][i];

        // add folded neutron charge distribution
        chd_tot[i] += chdf_np_vec[0][i];

        // add spin-orbit correction
        chd_ls[i] *= - std::pow( hbar_over_mc, 2 );
        chd_tot[i] += chd_ls[i];
    }

    // Find normalization of Charge Density Distribution
    double chd_norm = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        chd_norm += 4 * M_PI * chd_tot[i] * std::pow( rmesh[i], 2 ) * rdelt;
        
    }
    std::cout << "chd_norm = " << chd_norm << std::endl;

    // Write Charge Density Information
    double norm_fac = Nu_p.Z / chd_norm;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        chd_file << rmesh[i] << " " << chd_tot[i] * norm_fac << " " 
                 << chdf_np_vec[1][i] * norm_fac << " " 
                 << chdf_np_vec[0][i] * norm_fac << " " 
                 << chd_ls[i] * norm_fac << " " 
                 << chd_ls_test[i] * norm_fac  << " " 
                 << coulomb_vec[i] << std::endl;
    }

    chd_file.close();
    chd_file.clear();

    // Calculate RMS Radius
    double sum_pf = 0; // contribution from protons
    double sum_pfr2 = 0;
    double sum_nfr2 = 0; // contribution from neutrons
    double sum_lsr2 = 0; // contribution from spin-orbit correction
    double sum_totr2 = 0; // total
    for ( unsigned int i = 1; i < rmesh.size() - 1; ++i ) {
        
        sum_pf += std::pow( rmesh[i], 2 ) * chdf_np_vec[1][i] * rdelt;
        sum_pfr2 += std::pow( rmesh[i], 4 ) * chdf_np_vec[1][i] * rdelt;

        sum_nfr2 += std::pow( rmesh[i], 4 ) * chdf_np_vec[0][i] * rdelt;

        sum_lsr2 += std::pow( rmesh[i], 4 ) * chd_ls[i] * rdelt;

        sum_totr2 += std::pow( rmesh[i], 4 ) * chd_tot[i] * rdelt;

    }

    double rms_pf = std::sqrt( sum_pfr2 / sum_pf );
    double rms_tot = std::sqrt( 4 * M_PI * sum_totr2 / chd_norm );
    double rms2_p = 4 * M_PI * sum_pfr2 / chd_norm;
    double rms2_n = 4 * M_PI * sum_nfr2 / chd_norm;
    double rms2_ls = 4 * M_PI * sum_lsr2 / chd_norm;
    double rms_tot_check = std::sqrt( rms2_p + rms2_n + rms2_ls );

    std::string rms_filename = output_dir + nucleus_string + "_rms.out";
    std::ofstream rms_file( rms_filename.c_str() );
    rms_file << "Folded with proton charge distribution only" << std::endl;
    rms_file << "RMS radius: " << rms_pf << std::endl;
    rms_file << "Number of protons: " << sum_pf * 4 * M_PI << std::endl;

    rms_file << " " << std::endl;

    rms_file << "Proton, neutron and spin-orbit contributions included" 
             << std::endl;
    rms_file << "RMS radius: " << rms_tot << std::endl;
    rms_file << "Number of protons: " << chd_norm << std::endl;

    rms_file << "Mean Square Radius for different contributions "
             << "(normalized by Z = " << chd_norm << ")" << std::endl;

    rms_file << " " << std::endl;
    rms_file << "Proton Contribution (with Folding): " << rms2_p << std::endl;
    rms_file << "Neutron Conribution: " << rms2_n << std::endl;
    rms_file << "Spin-orbit Contribution: " << rms2_ls << std::endl;
    rms_file << " " << std::endl;
    rms_file << "Check with total: " << rms_tot_check << std::endl;

    rms_file.close();
    rms_file.clear();

    return 1;
}
