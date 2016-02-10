/***************************
  Have to write in terms of
  inputs and output the dom potential
  as a matrix in k space ?
  Probably better if put in Eigen here
 ***************************/
#include "dom_2_kSpace.h"
using namespace std;
int main()
{    
   double J;
   int n=0;
   int l=1;
   J=.5;

//  vector<double> r(96);
   double Ecm;//=1.;//=.2;//= 17.1376;//9.6573;// 19.5072;//29.2571;//10.7302;//20.4823;// 178.243;//63.3621;
   cin>>Ecm;
// We are using Bob's reduced mass A/(A+1) * 931.5
   double Mu = 914.79;//40.0/(1.+40.0)*m0;
 //  cout<<"MU   "<<Mu<<endl; 
 //  cout<<"Ecm   "<<Ecm<<endl; 
//cout <<n<<" l= "<<l<<" J=  "<<J<<endl;
   double ko;
   int nkpoints;
   vector<double> k,dk,Go;
 
   kgrid kpoints(Ecm,Mu); 
   kpoints.makeKmesh(k,dk);
   nkpoints =k.size();
   
   /*  r mesh  */
   double const rmax = 12.;
   int const  rpts = 200;
   
   double rdelt = rmax/rpts;
   vector<double> rmesh, dr;
   dr.assign(200,rdelt);
   for(int i =0; i<rpts; ++i)
       rmesh.push_back((i + 0.5)*rdelt );
/*   double ri = .01;  
   double rf = 10;  
   int Npoints = 200;
   double deltar = (rf-ri)/Npoints;
   vector< double > rmesh;
   for ( int i = 0; i <Npoints+1 ; ++i ) {
       rmesh.push_back( .01+ i * deltar   ); 
    }
*/
//////////////////////////////////////
//S_matrix caslculation///////////////
//////////////////////////////////////
/*      complex<double> S_lj,cS_lj;
     
      S_lj =complex<double>(0.,0.);
      cS_lj =complex<double>(0.,0.);


//    for(int jj = 0;jj<k.size();++jj){cout<<k[jj]<<endl;}
      cout<<k[0]<<"\t "<<k[k.size()-1] <<endl;
      kpoints.getPropagator(k,Go);
      ko = kpoints.getPole();
      for (int ii = 0 ; ii<12; ++ii){
         for (int jj = -1; jj<2; jj=jj+2){
             J = ii+jj/2.; 
             if (J<0 )    continue;                                        
             ReducibleSelfEnergy nonlocalDOM(n,ii,J,Ecm,Mu,k,dk,Go);

            //S matrix///
            ////////////
            nonlocalDOM.SelfEnergyAtPole(S_lj,cS_lj);

            //Calculating phase shift/// 
            //// *************************
            double phaseshift = atan2(imag(S_lj),real(S_lj));
            cout<<ii<<"\t"<<J<<"\t"<<real(S_lj)<<"\t"<<imag(S_lj)<<
                "\t"<<pow(real(S_lj),2)+pow(imag(S_lj),2)<<"\t phase shift= "<< phaseshift/2.0 <<endl;

//          cout<<"Ecm = " << Ecm <<"\t"<<"L = "<<ii<<"\t"<<J<<"\t"<<"arg S = "<<arg(S_lj)/2.<<"\t"<<
//                "real S = "<<real(S_lj)<<"\t"<<"Imagin S = " <<imag(S_lj)<<"\n";
        }
     }
*/
   p_spectralFunction dom_nonlocal(Mu,n,l,J);
//   dom_nonlocal.get_ParticleSpectralFunction2();
   dom_nonlocal.savespectral(Ecm);

//   dom_nonlocal.Sp_rr(Ecm,rmesh);
  // dom_nonlocal.Sp_rr(Ecm,rmesh,dr);
    
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
   return 0;
}
