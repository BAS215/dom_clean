#include "clebsch_gordan.h"  
int main() {
///////////////////////////////////////////
 ////////////////////initiation/////
 int lmax;
 int indexx;
 double theta;
 int theta_points = 30;
 string L_array[8] = { "s", "p", "d", "f", "g", "h", "i", "j" };
 string J_array[8] = { "1", "3", "5", "7", "9", "11", "13", "15" };

 double m_s1 = .5;   
 double m_s2 = .5;  
 double spin = .5;

///////////////////////////////////////////
 /////////////////////////////kmesh
 //
 vector<double> kmesh;
 vector<double> kweights;
 double const kmax = 6.0;
 int const kpts = 104;
 kmesh.resize( kpts );
 kweights.resize( kpts );
 GausLeg( 0.01, kmax, kmesh, kweights );
///////////////////////////////////////////
////////////////////initiation of vectors/////
vector<matrix_t> nkmatrixx = NKK(kmesh); 
cout<<"JDJJDJDJDJD"<<nkmatrixx.size()<<endl;
  
///////////////////////////////////////////
///////writing the final result into the output file
///////////////////////////////////////////
//
 string output_dir = "Output_test/New/nofk/";
 string outfilename = output_dir +"nk_total.out";

 std::ofstream outfile( outfilename.c_str() );

 cout<<"writing to the output file " <<endl;
 cout<<"theta_point " << theta_points<<endl;

 //for (int th = 0;th < theta_points; th++){
      theta = M_PI/3.; 
   //  theta =th * M_PI/theta_points;

     outfile<<theta<<endl;
     
     matrix_t F = FfF(kmesh, nkmatrixx, m_s1, m_s2, theta );

     for (int k1 = 0 ; k1<kmesh.size();++k1 ){
         
           for (int k2 = 0 ; k2<kmesh.size();++k2){
                   
               outfile <<kmesh[k1]<<" "<< kmesh[k2] << " "
                       << (F(k1,k2)) << " " 
                       << (nkmatrixx[0](k1,k2)) << " "
                       << (nkmatrixx[1](k1,k2)) << " "
                       << (nkmatrixx[2](k1,k2)) << " "
                       << (nkmatrixx[4](k1,k2)) << " "
                       << (nkmatrixx[5](k1,k2)) << " "
                       << (nkmatrixx[6](k1,k2)) << " "
                       << (nkmatrixx[7](k1,k2)) << " "
                       << (nkmatrixx[8](k1,k2)) << " "
                       << (nkmatrixx[9](k1,k2)) << " "
                       << (nkmatrixx[10](k1,k2)) << endl;  

           }     
       }
//  }
outfile.close();
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
return 0;
}
