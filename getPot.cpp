//
#include "pot.h"
#include "read_parameters.h"
#include "types.h"


int main( ) {

    // read in tz
    std::cout << "Enter tz (-0.5 or +0.5): " << std::endl;

    double tz;
    std::cin >> tz;

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

    config_file.close();
    config_file.clear();

    std::string output_dir = "Output_" + parameters_string + "/Output_vint/";
    std::string parameters_filename = parameters_string + ".inp";

    std::cout << "rmax = " << rmax << std::endl;
    std::cout << "rpts = " << rpts << std::endl;
    std::cout << "lmax = " << lmax << std::endl;

    // Create radial grid
    std::vector<double> rmesh;
    std::vector<double> rweights;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        rweights.push_back( rdelt );
        
    }

    // Create energy grid
    double Emax = 150;
    double Emin = -150;
    double edelt = 2;
    int epts = static_cast<int>( ( Emax - Emin ) / edelt );

    std::vector<double> emesh;
    for ( int i = 0; i < epts; ++i ) {

        emesh.push_back( Emin + i * edelt );
    }

    std::string tz_str;
    if ( tz > 0 ) tz_str = "p" + nucleus_string;
    else tz_str = "n" + nucleus_string;

    std::string input_dir = "Input/";
    std::string tz_filename = input_dir + tz_str + ".inp";

    // Create Nuclear Parameter Objects
    NuclearParameters Nu = read_nucleus_parameters( tz_filename );

    // Construct Parameters Object
    double Zp = tz + 0.5; // 1 for protons, 0 for neutrons
    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

    //Potential specifications
    int type = 1; // 1 is P.B. form (average), 0 is V.N. form
    int mvolume = 4;
    int AsyVolume = 1;

    // Construct Potential Object
    pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );

    std::string filename1 = output_dir + "coulomb.out";
    std::ofstream file1( filename1.c_str() );

    std::string filename4 = output_dir + "nonlocal_pot.out";
    std::ofstream file4( filename4.c_str() );

    int test_l = 0;
    double test_j = 0.5;
    double test_energy = Nu.Ef - 20;
    U.setEnergy( test_energy );
    U.setAM( test_l, test_j );

    for( unsigned int i = 0; i < rmesh.size(); ++i ) {

       //file1 << rmesh[i] << " " << U.coulomb_homogenous( rmesh[i] )
       //                  << " " << U.coulomb_exp( rmesh[i] ) << std::endl;

       for ( unsigned int j = 0; j < rmesh.size(); ++j ) {

            file4 << rmesh[i] << " " << rmesh[j] << " " 
                  << U.nonlocal_surface_IM( rmesh[i], rmesh[j] )
                  << std::endl;
       }
    }

    file1.close();
    file1.clear();
    file4.close();
    file4.clear();

    std::string filename2 = output_dir + "imaginary.out";
    std::ofstream file2( filename2.c_str() );

    std::string filename3 = output_dir + "dispersive.out";
    std::ofstream file3( filename3.c_str() );

    for ( unsigned int n = 0; n < emesh.size(); ++n ) {
        
        std::complex< double > vol_pot = U.volumeE( emesh[n] );
        std::complex< double > surf_pot = U.surfaceE( emesh[n] );

        file2 << emesh[n] - Nu.Ef << " " << imag( vol_pot ) 
              << " " << imag( surf_pot ) << std::endl;

        file3 << emesh[n] - Nu.Ef << " " << real( vol_pot ) 
              << " " << real( surf_pot ) << std::endl;
    }

    file2.close();
    file2.clear();
    file3.close();
    file3.clear();

    // radial dependence of HF

    return 0;
}
