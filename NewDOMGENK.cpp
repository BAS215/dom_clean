//
/*#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <gsl/gsl_sf_bessel.h>
#include "read_parameters.h"
#include "legendre.h"
#include "boundRspace.h"
//#include "manybody.h"
#include "boundRspace.h"
#include "meshes.h"
//#include "numerical.h"
#include "io.h"
#include "Gauss_q.h"
*/

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


///////////////////////////

    // Loop over protons and neutrons
    double Zp; // number of protons in projectile
    for ( unsigned int nu = 0; nu < Nu_vec.size(); ++nu ) {
        
        double tz = nu - 0.5; // +0.5 for protons, -0.5 for neutrons
        double Elim;
        if ( tz < 0 ) {

            Zp = 0;
         //   Elim = -62; // below 0s1/2 level for neutrons in 40Ca
        }
        else { 
            Zp = 1;
        //    Elim = -56; // below 0s1/2 level for protons in 40Ca
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



        // Loop over lj channels
        for ( int l = 0; l < lmax + 1; ++l ) {

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

                } // end loop over energy

            } //end loop over j
        } // end loop over l
    } // end loop over tz
    return 1;
}
