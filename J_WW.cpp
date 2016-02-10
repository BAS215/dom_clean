#include"pot.h"
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
int index_from_LJ( int L, double J ) {
    int Lmax = 5;
    // check to make sure that L is not greater than Lmax
    if ( L > Lmax ) {

        std::cout << "in function 'index_from_LJ': "
                  << "L too large." << std::endl;

        std::cout << "L = " << L << std::endl;
        std::cout << "Lmax = " << Lmax << std::endl;
    }

    int index;
    if ( J > static_cast<double>( L ) ) index = 2 * L;
    else index = 2 * L - 1;

    return index;

}
int main()
{

    std::string input_dir = "Input/";

    // Read Configuration file
    std::string config_filename = "ca40.config";
    std::ifstream config_file( config_filename.c_str() );

    std::string parameters_string;
    int fit_ph;
    double rmax = 12.;
    int rpts;
    int lmax;

    config_file >> parameters_string >> fit_ph;
    config_file >> rmax >> rpts >> lmax;

    int num_lj;
    config_file >> num_lj;

    config_file.close();
    config_file.clear();

 std::string output_dir = "Output_test/J_W/";
  //  std::string output_dir = "Output_" + parameters_string + "/J_W/";
  //  std::string parameters_filename = parameters_string + ".inp";
    std::string parameters_filename = "hosfitJune5.inp";

    std::cout << "rmax = " << rmax << std::endl;
    std::cout << "rpts = " << rpts << std::endl;
    std::cout << "lmax = " << lmax << std::endl;

    std::string n_filename = "Input/nca40.inp";
    std::string p_filename = "Input/pca40.inp";

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


    //Potential specifications
    int type = 1; // 1 is P.B. form (average), 0 is V.N. form
    int mvolume = 4;
    int AsyVolume = 1;

    /* CALCULATIONS */

    double Zp = 1.;
    double tz = .5; 
    // Nucleus Object
    const NuclearParameters &Nu = Nu_vec[1];

   // Construct Parameters Object
    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

    // Construct Potential Object
    pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
    pot *U1 = &U;
 
    // Create radial grid
    std::vector<double> rmesh;
    double rdelt = rmax / rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
       
    }
    // Create Energy Grid e
    double Emin = -300.;
    double Emax = 300;
    double deltaE = 5.;
    int epts = static_cast<int>( ( Emax - Emin ) / deltaE ) + 1;
    std::vector<double> emesh;
//    for ( int i = 0; i < epts; ++i ) {

  //         emesh.push_back( Emin + i * deltaE );
    //    }
    
/*    std::string outii_totb_old = output_dir + "J_W_tot_below_old.out"; 
    std::string outii_tota_old = output_dir + "J_W_tot_above_old.out"; 
    std::ofstream out_tot_below_old( outii_totb_old.c_str() );
    std::ofstream out_tot_above_old( outii_tota_old.c_str() );

    std::string outii_Below_old = output_dir + "J_W_Below_old.out"; 
    std::string outii_Above_old = output_dir + "J_W_Above_old.out"; 

    std::ofstream outputfile_Below_old( outii_Below_old.c_str() );
    std::ofstream outputfile_Above_old( outii_Above_old.c_str() );


    std::string outiis_Below_old = output_dir + "J_Ws_Below_old.out"; 
    std::string outiis_Above_old = output_dir + "J_Ws_Above_old.out"; 

    std::ofstream outputfiles_Below_old( outiis_Below_old.c_str() );
    std::ofstream outputfiles_Above_old( outiis_Above_old.c_str() );

    std::string outEfile = output_dir + "outEfile.out"; 
    std::ofstream potEfile( outEfile.c_str() );

    std::string outEfile_sur = output_dir + "outEfile_sur.out"; 
    std::ofstream potEfile_sur ( outEfile_sur.c_str() );

    double sur_B_old [6];
    double sur_A_old [6];
    double vol_B_old [6];
    double vol_A_old [6];
*/
   // std::cout <<" " << U1->nonlocalHF(1,1)<<std::endl; 
    std::string vreal1 = output_dir + "vHF_3d.out"; 
    std::ofstream vreal1_out(vreal1.c_str() );

    std::string vreal2 = output_dir + "vHF_diagonal.out"; 
    std::ofstream vreal2_out( vreal2.c_str() );

    std::string vreal3 = output_dir + "vHF_r_intg.out"; 
    std::ofstream vreal3_out( vreal3.c_str() );

   matrix_t Jafar(rmesh.size(),rmesh.size());
   std::vector <matrix_t> U_lj;
   std::vector<double> U_1;
   std::vector<double> U_2;
   std::vector<double> U_3;
   std::vector<double> U_4;
   std::vector<double> U_5;
   std::vector<double> U_6;
   std::vector<double> U_7;
   std::vector<double> U_8;
   std::vector<double> U_9;
   std::vector<double> U_10;
   U_1.assign( rmesh.size(), 0 );
   U_2.assign( rmesh.size(), 0 );
   U_3.assign( rmesh.size(), 0 );
   U_4.assign( rmesh.size(), 0 );
   U_5.assign( rmesh.size(), 0 );
   U_6.assign( rmesh.size(), 0 );
   U_7.assign( rmesh.size(), 0 );
   U_8.assign( rmesh.size(), 0 );
   U_9.assign( rmesh.size(), 0 );
   U_10.assign( rmesh.size(), 0 );
   int index;
 //  for ( int l = 0; l <  1; ++l ) {
 //      for( int up = -1; up < 2; up+=2 ) {

 //           double xj = l + up / 2.0;
 //           if ( xj < 0 ) continue;
 //           index = index_from_LJ( l, xj );
 //           U1->setAM(l,xj);
            U1->setEnergy(100.);          
            int kir = 0;
            U1->setAM(kir,.5);
            for (int k= 0 ; k< rmesh.size(); ++k) {
                 for (int kk= 0 ; kk< rmesh.size(); ++kk) {
	        if(k == kk) { vreal2_out << rmesh[k]<<" "<<  U1->nonlocalHF(rmesh[k],rmesh[kk])<<std::endl;} 
         	 //   Jafar(k,kk) =  U1->nonlocalHF(rmesh[k],rmesh[kk]); 
                 }
            }
            U_lj.push_back(Jafar);
 //      }
 //  }
	std::cout<<index<<std::endl;           
   for (int k= 0 ; k< rmesh.size(); ++k) {
        for (int kk= 0 ; kk< rmesh.size(); ++kk) {

                U_1[k] += U_lj[0](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
            /*    U_2[k] += U_lj[2](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_3[k] += U_lj[3](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_4[k] += U_lj[4](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_5[k] += U_lj[5](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_6[k] += U_lj[6](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_7[k] += U_lj[7](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_8[k] += U_lj[8](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_9[k] += U_lj[9](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
                U_10[k] += U_lj[10](k,kk) * rmesh[kk] * rmesh[kk] * rdelt;
               */
////diagonal vHF(r,r)
//	        if(k == kk) { vreal2_out << rmesh[k]<<" "<< U_lj[0](k,kk)<<std::endl;} 
/*                        <<  U_lj[1](k,kk)  << " " <<  U_lj[2](k,kk)  << " " <<  U_lj[3](k,kk)  << " " <<  U_lj[4](k,kk)<<" "    
                        <<  U_lj[5](k,kk)  << " " <<  U_lj[6](k,kk)  << " " <<  U_lj[7](k,kk)  << " " <<  U_lj[8](k,kk)<<" " 
		        <<  U_lj[9](k,kk)  << " " <<  U_lj[10](k,kk)<< std::endl;} 
*/
	       vreal1_out  << rmesh[k]<<" "<<rmesh[kk]<< " "<<U_lj[0](k,kk)<<std::endl; 
/*                        <<  U_lj[1](k,kk)  << " " <<  U_lj[2](k,kk)  << " " <<  U_lj[3](k,kk)  << " " <<  U_lj[4](k,kk)<<" "    
                        <<  U_lj[5](k,kk)  << " " <<  U_lj[6](k,kk)  << " " <<  U_lj[7](k,kk)  << " " <<  U_lj[8](k,kk)<<" " 
		        <<  U_lj[9](k,kk)  << " " <<  U_lj[10](k,kk)<< std::endl; 
*/
          }
          vreal3_out <<rmesh[k]<<" "<<U_1[k]<<std::endl;
		          //     " "<<U_2[k]<<" "<<U_3[k]
                          //     <<" "<<U_4[k]<<" "<<U_5[k]<<" "<<U_6[k]
                          //    <<" "<<U_7[k]<<" "<<U_8[k]<<" "<<U_9[k]<<" "<<U_10[k]<<std::endl;
    }
/*
    for (int kk = 0 ; kk < epts ; ++kk) {

	double vol_below_integral_old = 0.;
	double vol_above_integral_old = 0.;
	double sur_below_integral_old = 0.;
	double sur_above_integral_old = 0.;

        emesh.push_back( Emin + kk * deltaE );

        U1->setEnergy(emesh[kk]);
        potEfile << emesh[kk] << " " << imag(U1->volumeE( emesh[kk] )) <<" " << -emesh[kk] << " " << imag(U1->volumeE(-emesh[kk] )) << std::endl;
        potEfile_sur << emesh[kk] << " " << imag(U1->surfaceE( emesh[kk] )) <<" " << -emesh[kk] << " " << imag(U1->surfaceE(-emesh[kk] )) << std::endl;

        for ( int jj = 0; jj < lmax+1 ; ++jj) {


             sur_B_old[jj] = U1->Old_Im_pot_Integral(emesh[kk],rpts,rmesh,rdelt,jj,0)/Nu.A;
	     sur_below_integral_old += sur_B_old[jj];
///////////////

             sur_A_old[jj] =  U1->Old_Im_pot_Integral(-emesh[kk],rpts,rmesh,rdelt,jj,2)/Nu.A;
	     sur_above_integral_old += sur_A_old[jj];
///////////////


             vol_B_old[jj] =  U1->Old_Im_pot_Integral(emesh[kk],rpts,rmesh,rdelt,jj,1)/Nu.A;
	     vol_below_integral_old += vol_B_old[jj];
////////////// 


             vol_A_old[jj] =  U1->Old_Im_pot_Integral(-emesh[kk],rpts,rmesh,rdelt,jj,3)/Nu.A;
	     vol_above_integral_old += vol_A_old[jj];
	} 

        outputfiles_Above_old << -emesh[kk] << " " << -sur_A_old[0] << " " << -sur_A_old[1] << " " << -sur_A_old[2] << " "
			  << -sur_A_old[3] << " " << -sur_A_old[4] << " "<< -sur_A_old[5] <<" " << -sur_above_integral_old << std::endl;

        outputfiles_Below_old << emesh[kk] << " " << -sur_B_old[0] << " " << -sur_B_old[1] << " " << -sur_B_old[2] << " "
			  << -sur_B_old[3] << " " << -sur_B_old[4] << " "<< -sur_B_old[5] <<" " << -sur_below_integral_old << std::endl;
////////////// 


        outputfile_Below_old << emesh[kk] << " " << -vol_B_old[0] << " " << -vol_B_old[1] << " " << -vol_B_old[2] << " "
			  << -vol_B_old[3] << " " << -vol_B_old[4] << " "<< -vol_B_old[5] <<" " << -vol_below_integral_old << std::endl;

////////////// 

        outputfile_Above_old << -emesh[kk] << " " << -vol_A_old[0] << " " << -vol_A_old[1] << " " << -vol_A_old[2] << " "
			  << -vol_A_old[3] << " " << -vol_A_old[4] << " "<< -vol_A_old[5] <<" " << -vol_above_integral_old << std::endl;

        out_tot_below_old <<emesh[kk]<<" "<<" "<<-vol_below_integral_old-sur_below_integral_old<<std::endl;
        out_tot_above_old <<-emesh[kk]<<" "<<" "<<-vol_above_integral_old-sur_above_integral_old<<std::endl;
   }
outputfiles_Below_old.close();
outputfile_Below_old.close();
outputfiles_Above_old.close();
outputfile_Above_old.close();
out_tot_below_old.close();
out_tot_above_old.close();
potEfile_sur.close(); 
potEfile.close(); 
*/
vreal1_out.close();
vreal2_out.close();
vreal3_out.close();
return 0;
}
