// Contatins functions for reading in the various data
// ( right now there is just one function )

#include "read_data.h"

data_sets_t
read_volume_integral_data( const std::string &data_dir, 
                           const std::string &prefix,
                           const std::string &suffix ) {

    int Lmax = 4;

    data_sets_t set_vector;
    for ( int L = 0; L < Lmax + 1; ++L ) {

        std::string L_string = util::IntToStr( L );
        std::string filename = data_dir + prefix + L_string + suffix;

        std::ifstream file( filename.c_str() );

        if ( file.is_open() ) {
            
            std::cout << filename << " opened" << std::endl;
        }
        else {

            std::cout << "could not open file " << filename << std::endl;
            throw std::exception();
        }

        int num_pts;
        file >> num_pts;

        std::vector<double> energies;
        std::vector<double> values;
        std::vector<double> errors;

        for ( int i = 0; i < num_pts; ++i ) {

            double energy;
            double value;
            double error;

            file >> energy >> value >> error; 

            energies.push_back( energy );
            values.push_back( value );
            errors.push_back( error );

        } // end loop over i

        file.close();
        file.clear();

        Data_set D( energies, values, errors );

        set_vector.push_back( D );
    } // end loop over L

    return set_vector;
}

// Reads in \int{dr'r'^2 Im\Sigma(r, r')}, for example
Data_set
read_general_set( const std::string filename ) {

    std::list< std::string > list1 = util::read_commented_file( filename );

    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> errors;

    BOOST_FOREACH( std::string line, list1 ) {

        std::vector< std::string > vec = util::split( line );

        double x = boost::lexical_cast< double >( vec.at( 0 ) );
        double y = boost::lexical_cast< double >( vec.at( 1 ) );
        double err = boost::lexical_cast< double >( vec.at( 2 ) );

        x_values.push_back( x );
        y_values.push_back( y );
        errors.push_back( err );
    }

    Data_set D( x_values, y_values, errors );

    return D;
}

Data_set2
read_matrix ( const std::string filename ) {

    std::ifstream file( filename.c_str() );
    if ( file.is_open() != 1 ) {

        std::cout << "cannot open file " << filename << std::endl;
        throw std::exception();
    }

    // Read in matrix size (assuming square matrix)
    int mtx_size;
    double error;
    file >> mtx_size;
    file >> error;

    std::cout << mtx_size << std::endl;
    std::cout << error << std::endl;

    // Read in mesh
    std::vector<double> mesh, error_vec;
    for ( int i = 0; i < mtx_size; ++i ) {

        double k;
        file >> k;
        mesh.push_back( k );
        error_vec.push_back( error );
    }

    // Read matrix values
    matrix_t mtx( mtx_size, mtx_size );
    for ( int i = 0; i < mtx_size; ++i ) {
    for ( int j = 0; j < mtx_size; ++j ) {

        double y;
        file >> y;
        mtx( i, j ) = y;

    }
    }

    return Data_set2( mesh, mtx, error_vec );
}

Data_set2
read_matrix_ratio ( const std::string filename1, const std::string filename2 ) {

    Data_set2 set1 = read_matrix( filename1 );
    Data_set2 set2 = read_matrix( filename2 );

    // make sure the meshes are the same
    const std::vector< double > &mesh1 = set1.mesh;
    const std::vector< double > &mesh2 = set2.mesh;

    if( mesh1.size() != mesh2.size() ) {
        std::cout << "In read_matrix_ratio: Meshes not the same size." 
                  << std::endl;
        throw std::exception();
    }

    for ( unsigned int i = 0; i < mesh1.size(); ++i ) {

        if ( mesh1[i] != mesh2[i] ) {
            
            std::cout << "In read_matrix_ratio: Meshes not the same values." 
                      << std::endl;
            throw std::exception();
        }
    }

    // construct ratio of matrices
    matrix_t mtx_ratio( mesh1.size(), mesh1.size() );
    for ( unsigned int i = 0; i < mesh1.size(); ++i ) {
    for ( unsigned int j = 0; j < mesh1.size(); ++j ) {

        if ( set1.mtx( i, j ) == 0 ) {
            std::cout << "dividing by zero." << std::endl;
            throw std::exception();
        }
        mtx_ratio( i, j ) = set2.mtx( i, j ) / set1.mtx( i, j ); 
    }
    }

    return Data_set2( mesh1, mtx_ratio, set1.errors );
}

Data_vec 
read_all_data( std::string reactions_filename ) {
    // Read in the reactions that will be used in the fit
    // and create data objects
    std::ifstream reactions_file( reactions_filename.c_str(), std::ios::in );

    int num_reactions;
    reactions_file >> num_reactions; // number of reactions being fitted

    std::cout << "num_reactions = " << num_reactions << std::endl;
    // create vector to hold all the Data objects

    std::string reaction_title;

    // Create a vector to hold a data struct for each reaction
    Data_vec All_Data;
    for ( int i = 0; i < num_reactions; i++ ) {
        
        reactions_file >> reaction_title;

        std::string reaction_filename = reaction_title + ".inp";
        std::ifstream reaction_file( reaction_filename.c_str() );

        // Read in nucleus information
        double Zp;
        double Z;
        double A;
        double Ef;
        reaction_file >> Zp >> Z >> A >> Ef;

        reaction_file.close();
        reaction_file.clear();

        double tz;
        if (Zp == 0 ) tz = -0.5;
        else tz = 0.5;


        // --- READ DATA --- //
        
        //open data file
        std::string data_dir = "Output/N7_hw10_FULL/";

        //___________________________________________________________________
        //  read in CDBonn volume integrals 
        //___________________________________________________________________

        // 2h1p contributions
/*
        // k-space data
        std::string jw_h_prefix = reaction_title + "_jw_";
        std::string jw_h_suffix = "_hole.dat";

        std::string jd_h_prefix = reaction_title + "_jd_sub_";
        std::string jd_h_suffix = "_2h1p.dat";
        // r-space data
        std::string jw_h_prefix = reaction_title + "_jw_";
        std::string jw_h_suffix = "_2h1p_r_space.dat";

        std::string jd_h_prefix = reaction_title + "_jd_sub_";
        std::string jd_h_suffix = "_2h1p_r_space.dat";

        data_sets_t jw_h_sets =
            read_volume_integral_data( data_dir, jw_h_prefix, jw_h_suffix );

        data_sets_t jd_h_sets =
            read_volume_integral_data( data_dir, jd_h_prefix, jd_h_suffix );

*/
        // 2p1h contributions
/*
        // k-space data
        std::string jw_p_prefix = reaction_title + "_jw_";
        std::string jw_p_suffix = "_particle.dat";

        std::string jd_p_prefix = reaction_title + "_jd_sub_";
        std::string jd_p_suffix = "_2p1h.dat";
        // r-space data
        std::string jw_p_prefix = reaction_title + "_jw_";
        std::string jw_p_suffix = "_2p1h_r_space_lowE.dat";

        std::string jd_p_prefix = reaction_title + "_jd_sub_";
        std::string jd_p_suffix = "_2p1h_r_space.dat";

        data_sets_t jw_p_sets =
            read_volume_integral_data( data_dir, jw_p_prefix, jw_p_suffix );

        data_sets_t jd_p_sets =
            read_volume_integral_data( data_dir, jd_p_prefix, jd_p_suffix );

*/

        // Integrated Form
        std::string int_form_filename1 = 
            data_dir + "re_shape_AV18_nca40_s_at_Ef_int.dat";

        Data_set Int_form_set = read_general_set( int_form_filename1 );

        // Diagonal Form
        std::string diag_form_filename1 =
            data_dir + "re_shape_AV18_nca40_s_at_Ef_diag.dat";

        Data_set Diag_form_set = read_general_set( diag_form_filename1 );

        // Create Data Struct
        Data D( A, Z, tz, Ef, Int_form_set, Diag_form_set );

        All_Data.push_back( D );

    } // end loop over reactions

    return All_Data;
}

