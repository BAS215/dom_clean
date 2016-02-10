#include "reaction.h"
#include "pot.h"
#include "boundRspace.h"
//#include <gsl/gsl_sf_legendre.h>
//#include <gsl/gsl_complex.h>
#include "read_parameters.h"
#include <string>
#include <sstream>
#include <fstream>

std::string
get_lj_string( const std::string L_string, double J ) {

    // the addition of 0.5 is to round to the nearest integer
    // instead of truncating the entire decimal part
    int J_int = static_cast< int >( 2 * J + 0.5 );

    std::ostringstream lj_strm;
    lj_strm << L_string << J_int;

    return lj_strm.str();
}

int main() {

  //reading boundrspace data to calculate chi squared//
  chd_data Chisq_data_exp;
  std::ifstream filein("ca40_bound_chd.dat");
  filein >> Chisq_data_exp.r_0 >> Chisq_data_exp.ch_den_0_exp
         >> Chisq_data_exp.err_0_exp >> Chisq_data_exp.r_middle_exp
         >> Chisq_data_exp.ch_den_middle_exp >> Chisq_data_exp.err_middle_exp
         >> Chisq_data_exp.R_rms_exp >> Chisq_data_exp.err_R_rms_exp;

    filein.close();

    // this is used for creating the names of the output files
    std::string L_array[8] = { "s", "p", "d", "f", "g", "h", "i", "j" };

//   double Zp = 0;
   double Zp = 1;
    int Lmax = 5;
   // int Lmax = 6;
    double rmax = 12.0;
    int rpts = 180;

    std::string np_str = "pca40";
   // std::string np_str = "nca40";
    std::string input_file = "Input/" + np_str + ".inp";
    std::string output_dir = "Output_test/Input/";
   // std::string parameters_filename ="InputCa40.inp";
    std::string parameters_filename ="Input19Feb2015.inp";
   // std::string parameters_filename ="Input.inp";

    // Read parameters
    NuclearParameters Nu = read_nucleus_parameters( input_file );
    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );


    // construct Potential Object
    int typeN = 1;
    int mvolume = 4;
    int AsyVolume = 1;
    double tz = 0.5;
    if ( Zp < 1 ) tz = -0.5;


    pot U0 = get_bobs_pot2( typeN, mvolume, AsyVolume, tz, Nu, p );
    pot *U1 = &U0;

    double Emax;

    // Construct boundRspace Object
     //Emax = 2. * E.;//  - Nu.Wgap * p.fGap_B; // energy < Ef where imaginary part begins
    //if (tz>0) Emax = -12. ; else 
  //  Emax =  2 * Nu.Ef;
    double Emin = -200 + Emax;
    Emax = 2 * Nu.Ef;
    boundRspace B( rmax, rpts, Nu.Ef, Nu.ph_gap, Lmax, Nu.Z, Zp, Nu.A, U1 );

		
    // Create radial grid
    std::vector< double > rmesh= B.make_rmesh();
    std::vector< double > rmesh_pN= B.make_rmesh();
    std::vector< double > delta_pN;//= B.make_rmesh();

    // Get bound levels
    double Efermi=Nu.Ef;
    double tol = 0.01;
    std::vector< lj_eigen_t > bound_levels = B.get_bound_levels( rmesh, tol );


    // Create energy meshes
    std::vector< mesh_t > emesh_vec = 
        B.get_emeshes( rmesh, Emin, Emax, bound_levels );


    // Calculate propagators
    std::vector< prop_t > prop_vec = B.get_propagators( rmesh, emesh_vec );


    // Calculate bound states and write them out to a file
    std::string J_array[8] = { "1", "3", "5", "7", "9", "11", "13", "15" };

    std::string qp_filename = 
            output_dir + np_str + "_quasiparticle" + ".out";
    std::ofstream qp_file( qp_filename.c_str() );


    for ( int L = 0; L < Lmax + 1; ++L ) {

        for ( int up = -1; up < 2; up+=2 ) {

            double J = L + up / 2.0;
            if ( J < 0 ) continue;

            int j_index = ( up - 1 ) / 2;
            std::string j_string;
            if ( L == 0 ) j_string = J_array[ L ];
            else j_string = J_array[ L + j_index ];

            std::string lj_string = get_lj_string( L_array[L], J ); 
        //    std::string lj_string = L_array[L] + j_string + "2";

            for ( int N = 0; N < B.max_n; ++N ) {

                double Estart = Nu.Ef;
                tol = 0.01;
                eigen_t bound_info = 
                    B.find_boundstate( rmesh, Estart, N, L, J, tol );

                double QPE = bound_info.first;
                std::vector< double > &QPF = bound_info.second;

                double S = B.sfactor( rmesh, QPE, L, J, QPF );


                if (QPE < 0 ){ qp_file << N << lj_string << " " << QPE 
                        << " " << S << " "<< std::endl;}

                if ( QPE > 0 ) {  

                  //  std::cout << "No bound state found for N = " 
                    //          << N << std::endl;
                    break; 
                }                        

            } // end loop over N


            // Calculate spectral function
            
            std::vector< double > s_of_E = B.spectral_strength( L, J, rmesh, emesh_vec, prop_vec ); 

 	    std::string spf_filename = output_dir + np_str + "_spectral_f_" + lj_string + "2.out";

            std::ofstream spf_file(spf_filename.c_str() );
    
   spf_file.close();
        } // end loop over J
    cout<<"end loop over J"<<endl;
    } // end loop over L



    // Calculate charge density
    std::vector< double > point_dist = 
        B.point_distribution( rmesh, Emax, emesh_vec, prop_vec, bound_levels );

    std::vector< double > point_dist_normalized ; 


    // Find normalization of point distribution
    double point_norm = 0;
    double point_norm_pN = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        cout<<rmesh[i]<<"\t"<< U1->coulomb_homogenous(rmesh[i])<<"\t"<<U1->coulomb_exp(rmesh[i])<<endl;
        point_dist_normalized.push_back( point_dist[i] );
        point_norm += 4 * M_PI * point_dist[i] * std::pow( rmesh[i], 2 ) * B.rdelt; 
      /*  if (rmesh[i] < 4.8) {
                             rmesh_pN.push_back(rmesh[i]); 
                             delta_pN.push_back( B.rdelt ); 
                            }

        else               {
                            // point_norm += 4 * M_PI * point_dist[i] * std::pow( rmesh[i], 2 ) * B.rdelt * std::pow((39./40.),3); 
                             point_norm += 4 * M_PI * point_dist[i] * std::pow( rmesh[i], 2 ) * B.rdelt ; 
                             rmesh_pN.push_back( (39./40.) * rmesh[i] ); 
                             if (rmesh[i-1] > rmesh_pN[i]) { rmesh_pN[i-1] = rmesh_pN[i] ; } 
                             delta_pN.push_back( (39./40.) * B.rdelt ); 
                             cout << rmesh_pN[i] <<" " <<point_norm<<endl; 
                            }  
     */
    }
    double norm_fac = Nu.Z / point_norm;

////////////////// ////////////// ////////////// //////////////
    cout<<point_norm<<" "<<point_norm_pN <<"B.rdelt    "<<B.rdelt<<"    norm fac   "<< norm_fac<<endl;
///////////////// ////////////// ////////////// /////////
//
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    
        point_dist_normalized[i] *= norm_fac;
    }
    std::cout << "Finishded compare3d" << std::endl;
//**************************************************************************
    std::cout << "Folding density distribution" << std::endl;

    // Fold point charge density with neutron and proton charge distribution
    // FIXME can only use the same point distribution for neutrons and protons
    // if the nucleus has N = Z. Eventually, the code below needs to be made
    // more general.

    std::vector< double > rweights;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        rweights.push_back( B.rdelt );
        delta_pN.push_back( B.rdelt );
    }

    std::vector<double> proton_folding = 
        folded_ch_density( rmesh, rweights, point_dist, tz, Nu.A );

//    std::vector<double> proton_folding = 
//        folded_ch_density( rmesh_pN, delta_pN, point_dist, tz, Nu.A );


    std::vector< double > chdf; // folded charge density
    double chdfolded_norm=0.;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

     //   chdf.push_back( proton_folding[i] + neutron_folding[i] );
        chdf.push_back( proton_folding[i]  );
        chdfolded_norm += 4 * M_PI * chdf[i] * std::pow( rmesh_pN[i], 2 ) * delta_pN[i];
    }
    
    std::string chd_filename = output_dir + np_str + "_folded_unfolded_charge_density.out";
    std::ofstream chd_file( chd_filename.c_str() );
    
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        chd_file << rmesh[i]<<" "<<rmesh_pN[i] << " " << chdf[i] <<" "<<
              " "<< point_dist[i] << " "<< point_dist_normalized[i] << std::endl;
    }

    chd_file.close();
    chd_file.clear();



    // Miscellany
    std::string misc_filename = output_dir + np_str + "_misc.out";
    std::ofstream misc_file( misc_filename.c_str() );
   
 
   // Calculate root-mean-square radius

    double rms1 = B.chd_rms_radius( rmesh, chdf );
    double rms11 = B.chd_rms_radius( rmesh_pN, chdf );
    double rms = B.chd_rms_radius( rmesh, point_dist );
    double rms2 = B.chd_rms_radius( rmesh_pN, point_dist );

    cout<<rms1<<" "<<rms11<<" "<<rms<<" "<<rms2<<endl;

    misc_file << "rms radius:(chdf, point_dist) " <<rms1 << " "<<rms<<"\n" 
              << "rms experimental value :" << Chisq_data_exp.R_rms_exp << " "
              <<"+-"<<" "<<  Chisq_data_exp.err_R_rms_exp << std::endl;

    // Calculate Particle Number 
    misc_file << "Particle Number (point_dist normalization: " << point_norm <<"pN "<<point_norm_pN<<"\n"<< 
               "Particle number (chdf folded charge distribution normalization)"<<chdfolded_norm<<"\n"<<
               "Particle experimental value :" << 20 << std::endl;

    misc_file.close();
    misc_file.clear();

    return 0;

}
