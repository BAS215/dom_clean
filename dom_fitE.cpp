// This program is used when fitting the nonlocal parameters
// to the get the energy levels. It only calculates the
// quasiparticle and quasihole energies 
#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include "read_parameters.h"
#include "boundRspace.h"
#include "meshes.h"
#include "numerical.h"
#include "io.h"


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

    // Create radial grid
    std::vector<double> rmesh;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        
    }

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
    for ( unsigned int nu = 1; nu < 2; ++nu ) {
    //for ( unsigned int nu = 0; nu < Nu_vec.size(); ++nu ) {
        
        double tz = nu - 0.5; // +0.5 for protons, -0.5 for neutrons
        std::map< std::string, int > lj_map = map_vec[nu];
        if ( tz < 0 ) Zp = 0;
        else Zp = 1;

        // Nucleus Object
        const NuclearParameters &Nu = Nu_vec[nu];

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

        // Construct Potential Object
        pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
        pot *U1 = &U;

        // Construct boundRspace Object
        boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, lmax, Nu.Z, Zp, Nu.A, U1 );

        double tol = 0.01;
        std::vector< lj_eigen_t > bound_levels = 
            B.get_bound_levels( rmesh, tol );

        // Output Files
        std::string qp_filename = output_dir + np_strings[nu] 
                                + "_quasiparticle" + ".test";
        std::ofstream qp_file( qp_filename.c_str() );

        // Loop over lj channels
        for ( int l = 0; l < lmax + 1; ++l ) {

            for( int up = -1; up < 2; up+=2 ) {

                double xj = l + up / 2.0;
                int j_index = ( up - 1 ) / 2;
                if ( xj < 0 ) continue;

                int index = B.index_from_LJ( l, xj );
                std::vector< eigen_t > &bound_info = bound_levels[index];

                double lp; 
                if ( xj > l ) lp = l;
                else lp = - l - 1;
                
                std::string j_string;
                if ( l == 0 ) j_string = J_array[ l ];
                else j_string = J_array[ l + j_index ];

                std::string lj_string = L_array[l] + j_string + "2";

                /* QUASIPARTICLE AND QUASIHOLE INFORMATION */

                // Find self-consistent solutions

                // loop over the orbitals and write out the 
                // bound state information to a file
                int intj = static_cast<int>( 2 * xj );
                for ( unsigned int N = 0; N < bound_info.size(); ++N ) {

                    std::string level_str = util::IntToStr( N ) + L_array[l] 
                                          + j_string + "2";

                    std::string wfx_filename = output_dir + np_strings[nu] 
                        + "_" + level_str + ".wfx";
                    std::ofstream wfx_file( wfx_filename.c_str() );

                    // Quasiparticle energy
                    double QPE = bound_info[N].first; 

                    // Quasiparticle wavefunction
                    std::vector<double> &QPF = bound_info[N].second;
                    double S = B.sfactor( rmesh, QPE, l, xj, QPF );

                    qp_file << N << L_array[l] << intj << "/2" 
                            << " " << QPE << " " << S << std::endl;

                    // Write out wavefunctions to a file
                    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
                        
                        wfx_file << rmesh[i] << " " << QPF[i] << std::endl;

                    }

                    wfx_file.close();
                    wfx_file.clear();

                } // end loop over N 

            } //end loop over j
        } // end loop over l

        qp_file.close();
        qp_file.clear();

    } // end loop over Nu_vec

    return 1;
}
