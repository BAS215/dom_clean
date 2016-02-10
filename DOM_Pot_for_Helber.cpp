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

    std::string output_dir = "Output_test/DOM_pot_for_Helber/";
    std::string parameters_filename = "InputCa48.inp";

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
 
    //bound state properties for calculating  wavefunction
    boundRspace initiate(rmax ,180  , -4.7, 0 , 5 , 20 , Zp , 40 , U1);
  


    // Create momentum space grid
    std::vector<double> kmesh;
    std::vector<double> kweights;
    double const kmax = 6.0;
    int const kpts = 104;
    kmesh.resize( kpts );
    kweights.resize( kpts );
    GausLeg( 0., kmax, kmesh, kweights );

    // Create radial grid
    std::string r_grids = output_dir + "rmesh.out"; 
    std::ofstream rmesh_out( r_grids.c_str() );
    std::vector<double> rmesh;
    std::vector<double> rweights;
    double rdelt = rmax / rpts;

    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
        rweights.push_back( rdelt );

        rmesh_out<<rmesh[i]<<" "<<rweights[i]<<std::endl;        
    }

    rmesh_out.close();  

//**************************************************************************************
//    **************************************************
    //Non local potential for different energy values
    //**************************************************


    
  for (int en = 2 ; en<50 ;  en = en+20)
  {
    double En =en ;
    U1->setEnergy(-En);
    /////starting from L,J
    int index = 0; 

    int L = 2;
    int up = -1;		 
   // Create Bessel Function matrix in k and r
    matrix_t bess_mtx( kmesh.size(), rmesh.size() );
    for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {
        for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

        double rho = kmesh[nk] * rmesh[nr];
        bess_mtx( nk, nr ) = gsl_sf_bessel_jl( L, rho );
        }
    }

    cout << "L is "<< " " << L<<endl;
    double J = L + up / 2.0;
    U1->setAM(L,J);

    //converting "index" to string//

    //string result_index; 
    //ostringstream convert_index; 
    // convert_index << index;  
    //result_index = convert_index.str(); 
    //////////////////////////// 


    //Energy loop 0 to 100 step 2 for each L, J
    //

    //double en;          
    //cin>> en ;
    //converting "en " to string for labeling files appropriately//

    //string result_energy;   
    //ostringstream convert2; 
    //convert2 << en;      
    //result_energy = convert2.str(); 
    //////////////////////////// 


    cout<<"  "<< " energy file for E" << "  " <<En<< endl ; 
    //
    std::string nonlocalfile = "test50.out" ;// output_dir + "V_lj_index_" + result_index +  "_E_" + result_energy + ".out"; 
    std::ofstream nonE( nonlocalfile.c_str() );
            
    //std::string wave_file = output_dir + "wave_function" + result_energy + ".out"; 
    //std::ofstream wavef( wave_file.c_str() );
    //////////////////////////////////////////////////////////////////////////////////////////
    //wave function in r space
    //////////////////////////////////////////////////////////////////////////////////////////
    double tol=.05;
    double Ef = -4.7;
    double estart=Ef;
    //  boundRspace::find_boundstate( const std::vector<double> &rmesh, double Estart, 
    //                                int N, int l, double j, double tol ) {
    //
    eigen_t waves_function = initiate.find_boundstate(rmesh, estart,0,L,J,tol);

    //building V(r,r')and store it in V_r_matrix
    vector<double> WaveFunction ( rmesh.size() );
    cmatrix_t V_r_nonlocal_matrix( rmesh.size(), rmesh.size() );
    cvector_t V_r_local_matrix( rmesh.size());
    for( int i = 0; i < rpts; ++i ) {
        for( int j = 0; j < rpts; ++j ) {

            V_r_nonlocal_matrix(i,j) = complex<double>( 0.0, 0.0 );
            V_r_nonlocal_matrix(i,j) = U1->nonlocalPart( rmesh[i], rmesh[j]) ;

        }

       //building the local potential V(r)and store it in V_r_matrix
        V_r_local_matrix[i] = 0.0;
        WaveFunction[i]=0.0;
        V_r_local_matrix[i] = U1->localPart(rmesh[i]);
        WaveFunction[i]=waves_function.second[i];
               
    }  
    cout <<"for E = "<<" " << En<<" "<<"v_local at " <<rmesh[13]
        <<" "<<"is" <<" "<<V_r_local_matrix[13]<<" "<<U1->localPart(rmesh[13])<<endl;

    cout <<"for E = "<<" " << En<<" "<<"v_NONlocal at " <<rmesh[13]
         <<" "<<"is" <<" "<<V_r_nonlocal_matrix(12,13)<<" "<< U1->nonlocalPart( rmesh[12], rmesh[13])<< endl;
    //Fourier bessel transform of V(r,r') and save V(k,k') in the output file 
    //
    for( unsigned int nk = 1; nk < kmesh.size(); ++nk ) {

        complex<double> sum_local = 0.; 
        double sum_wave = 0.; 
       
        //bessel fourier transformation/integral of V(r) and wave function
        for( int ii = 0 ; ii <rmesh.size(); ++ii){ 
            sum_local = sum_local + bess_mtx( nk, ii ) * V_r_local_matrix[ii] * rmesh[ii] * rmesh [ii] * rdelt ;
            sum_wave = sum_wave + bess_mtx( nk, ii ) * WaveFunction[ii] * rmesh[ii] * rmesh [ii] * rdelt ;
        }

        //wavef<<kmesh[nk]<< " " <<sum_wave<<endl;                                                      

        for( unsigned int mk = 0; mk < kmesh.size(); ++mk ) {

            complex<double> sum_out = 0.; 

            for( int i = 0; i < rpts; ++i ) {

                double jl1 = bess_mtx( nk, i );
                complex<double> sum_in = 0.; 

                for( int j = 0; j < rpts; ++j ) {
                      
                    double jl2 = bess_mtx( mk, j );
                    sum_in = sum_in + V_r_nonlocal_matrix( rmesh[i], rmesh[j]) * jl2 * rmesh[j] * rmesh [j] * rdelt ;

                 }//inner r loop

                 sum_out = sum_out + sum_in * jl1 * rmesh[i] * rmesh[i] * rdelt;

             } //outer r loop

             double realPot = real(sum_out);  
             double imagPot =  imag(sum_out);  

             if ( nk == mk ) { 
                              realPot = realPot + real(sum_local); 
                              imagPot = imagPot + imag(sum_local); 
             }//if 

             nonE <<kmesh[nk]<<" "<<kmesh[mk]<< " "<< realPot << " " << imagPot <<std::endl;

        }//end loop over mk(second,inner,  loop over k)
   }//end loop over nk(first  loop over k)

   nonE.close();
   // wavef.close();
}//end of the energy loop
   return 0;
}
