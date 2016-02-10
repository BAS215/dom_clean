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

   std::string output_dir = "Output_test/J_W/";
  //  std::string output_dir = "Output_" + parameters_string + "/J_W/";
  //  std::string parameters_filename = parameters_string + ".inp";
    std::string parameters_filename = "InputCa40.inp";

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

    double Zp = 0;//1.;
    double tz = -.5;//.5; 
    // Nucleus Object
    const NuclearParameters &Nu = Nu_vec[0];

   // Construct Parameters Object
    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

    // Construct Potential Object
    pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
    pot *U1 = &U;
 
    //bound state properties for calculating  wavefunction
    boundRspace initiate(12 ,180  , -12.0, 0 , 5 , 20 , Zp , 40 , U1);
    // Create momentum space grid
    std::vector<double> kmesh;
    std::vector<double> kweights;
    double const kmax = 6.0;
    int const kpts = 104;
    kmesh.resize( kpts );
    kweights.resize( kpts );
    GausLeg( 0., kmax, kmesh, kweights );
    //************************************

    // Create radial grid
    std::vector<double> rmesh;
    double rdelt = rmax / rpts;
    for ( int i = 0; i < rpts; ++i ) {
       
       rmesh.push_back( ( i + 0.5 ) * rdelt ); 
    }
    //************************************
    
 
     matrix_t bess_mtx( kmesh.size(), rmesh.size() );
     cmatrix_t V_local_matrix(rpts,rpts);
     cmatrix_t V_nonlocal_matrix(rpts,rpts);
     cmatrix_t V_r_tot ;
     vector<double> wavefunction; 
///////////////////Non local potential for an Energy value on diagonal


    /////calculating rmesh and nonlocalPart(r,r) stored in the above output file
    
  //  for(int en = 2; en <4 ; en = en+2) { 
 //      int en; 
//       cin >> en;
 //      string result;   ostringstream convert; convert << en;   result = convert.str(); 
double En = -15.6;
//for (int ii=0;ii<3; ++ii){
      // En = En+15.;
       cout<<"E = "<< En<<endl;
       U1->setEnergy(En);
       U1->setAM(0,1.5);
       cout<<" SO V  "<< U1->SpinOrbit.V << endl; 
       cout<<" SO Vzero "<< U1->SpinOrbit.Vzero << endl; 
       cout<<" SO Vdisp "<< U1->SpinOrbit.V-U1->SpinOrbit.Vzero << endl; 
       cout<<" S_Above  "<< U1->SurfaceAbove.V << endl; 
       cout<<" S_Below  "<< U1->SurfaceBelow.V << endl; 
       cout<<" V_Above  "<< U1->VolumeAbove.Vvol << endl; 
       cout<<" V_Below  "<< U1->VolumeBelow.Vvol << endl; 
//}
       //////output file 
       std::string nonlocalfile = output_dir + "nonlocal_V_k_Ener_" + ".out"; 
       std::ofstream nonE( nonlocalfile.c_str() );

       std::string wfile = output_dir + "wavefunction_E_" + ".out"; 
       std::ofstream wavefile( wfile.c_str() );

       for(int L = 0; L < 1 ; L++) { 
            
           for ( int up = -1; up <= 1; up+=2 ) {
                double J = L + up / 2.0;
  
                if ( J < 0 ) continue;

                U1->setAM(L,J);
                               
                //preparing Potential matrices in r space
                for( int i = 0; i < rpts; ++i ) {
                     for( int j = 0; j < rpts; ++j ) {
                          V_nonlocal_matrix(i,j) = complex<double>( 0.0, 0.0 );
                     }
                     V_local_matrix(i,i) = complex<double>( 0.0, 0.0 );
                     wavefunction.push_back( 0.); 
                }

		//wavefunction
                eigen_t waves_function = initiate.find_boundstate(rmesh, -4.7,0,L,J,.05);
/*
               // Create Bessel Function matrix in k and r////////////"<<endl;
               for( unsigned int nk = 0; nk < kmesh.size(); ++nk ) {
                  for( unsigned int nr = 0; nr < rmesh.size(); ++nr ) {

   	             bess_mtx( nk, nr ) = 0.;
                     double rho = kmesh[nk] * rmesh[nr];
                     bess_mtx( nk, nr ) = gsl_sf_bessel_jl( L, rho );
                  }
                }
               
*/
                //*****************************************************
                

               //building the local and local potential in r space"<<endl;
               double Sum_nonlocal_real = 0.0;
               double Sum_nonlocal_imag = 0.0;
               for( int i = 0; i < rpts; ++i ) {
                     wavefunction[i] = waves_function.second[i]; 
                   //  cout<<rmesh[i]<<" "<< waves_function.second[i]<<endl; 
                     for( int j = 0; j < rpts; ++j)  {    

                         V_nonlocal_matrix(i,j) = (U1->nonlocalPart( rmesh[i], rmesh[j]));

			 //Integrating the real Part of the potentila over space                         
			 complex<double> Vij =  V_nonlocal_matrix(i,j);

                         Sum_nonlocal_real = Sum_nonlocal_real + rmesh[i] * rmesh[j] * real(Vij) * rdelt * rdelt; 
                         Sum_nonlocal_imag = Sum_nonlocal_imag + rmesh[i] * rmesh[j] * imag(Vij) * rdelt * rdelt; 

                         if (i==j)  V_nonlocal_matrix(i,j) += (U1->localPart( rmesh[i]));
                     }
               } 
               cout<<"L "<<L<<" J "<<J<< "\n Sum_nonlocal_real "<<Sum_nonlocal_real<< "\n Sum_nonlocal_imag "<<Sum_nonlocal_imag<<endl;   
               //*****************************************************
               //*****************************************************
               /*
               complex_t sum_out,sum_in; 


	       double sum_wave ;
               //Fourier transform of V(r,r') to V(k,k') 
               for(int m = 0 ; m<kmesh.size() ; ++m) {
                   for(int n = 0 ; n<kmesh.size() ; ++n) {
 		     sum_wave = 0.;
                     sum_out =  complex<double> (0.,0.);
                      for( int i = 0; i < rpts; ++i ) {

                         double jl1 = bess_mtx( m, i );
                         sum_in =  complex<double> (0.,0.);

                         for( int j = 0; j < rpts; ++j)  {    

                            double jl2 = bess_mtx( n, j );
                            sum_in = sum_in +  V_nonlocal_matrix(i,j) * jl2 * rmesh[j] * rmesh [j] * rdelt ;

                         }
                         sum_out +=  sum_in * jl1 * rmesh[i] * rmesh [i] * rdelt ;
                         sum_wave +=  wavefunction[i] * jl1 * rmesh[i] * rmesh [i] * rdelt; 
                      } 
                      nonE<<kmesh[m]<<" " <<kmesh[n] <<" " <<real(sum_out)<<" "<<imag(sum_out)<<endl;   
		      if (m == n) { wavefile<<kmesh[n]<<" "<<sum_wave<<endl;}
                   }//inner k, n
                }//outer k, m

*/
           }//end loop over up
       }//end loop over L

       nonE.close();
       wavefile.close();
//   }//end loop over Energy

return 0;

}
