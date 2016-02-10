/***************************
  Have to write in terms of
  inputs and output the dom potential
  as a matrix in k space ?
  Probably better if put in Eigen here
 ***************************/
#include"dom_2_kSpace.h"
using namespace std;
int main()
{    
   double Ecm, J;
   int n=0, l=5;
    J=5.5;

    vector<double> rmesh;
    double rmax = 12.;
    int rpts = 180;
    double rdelt = rmax/rpts;
    for( int i = 0; i < rpts; ++i ) {

        rmesh.push_back( ( i + 0.5 ) * rdelt );
    }
  // cout<<"n =";
  // cin>> n;
  // cout<<"l =";
  // cin>> l;
  // cout<<"J =";
  // cin>> J;

 //  int lineS0 = 141;
 //  double dat1,dat2,dat3;	
 //  vector<double> Ecm_vector;
 //  std::ifstream Ecmfile( "Output_hosfit/Smatrix/L4up.dat" );
 //  for (int i=0 ; i <lineS0 ; ++i){
 //       Ecmfile >> dat1 >> dat2 >> dat3 ;
 //       cout << dat1<<endl;
 //       Ecm_vector.push_back(dat1);
 //   }	 
  double dE;
  for (int ii = 0; ii<100; ++ii){
  Ecm = 1.;//10.73;
  //  /// TEsting what we have till now //////
 //   for (int ii=0; ii<Ecm_vector.size() ; ii++){
   // double Ecm = Ecm_vector[ii];
    double Mu = M;
    double ko;
    int nkpoints;
    vector<double> k,dk,Go;
    kgrid kpoints(Ecm,Mu); 
    kpoints.makeKmesh(k,dk);
    nkpoints =k.size();
    kpoints.getPropagator(k,Go);
    ko = kpoints.getPole();
    ////////////////////////////////////////////
    /*    On-shell calculation   */ 
  //   complex<double> S_lj,cS_lj;
     
  //     S_lj =complex<double>(0.,0.);
  //    cS_lj =complex<double>(0.,0.);
/////*    
////      cin>>l;
////      J = l + 0.5;
////*/    
       ReducibleSelfEnergy nonlocalDOM(n,l,J,Ecm,Mu,k,dk,Go);
       cout <<" MAss= " << Mu<<endl;
//       nonlocalDOM.SelfEnergyAtPole(S_lj,cS_lj);
//       cout<<"Ecm = " << Ecm <<"\t"<<"L = "<<l<<"\t"<<J<<"\t"<<"arg S = "<<arg(S_lj)/2.<<"\t"<<
//             "real S = "<<real(S_lj)<<"\t"<<"Imagin S = " <<imag(S_lj)<<"\n";
//  }

    /***** Full reducible self-energy ******/
/*
     l = 0;
     J = 0.5;
    vector<double> wf;
    MatrixXcd Sigma(nkpoints,nkpoints), cSigma(nkpoints,nkpoints);
*/
  //  ReducibleSelfEnergy nonlocalDOM(n,l,J,Ecm,Mu,k,dk,Go);
  //  nonlocalDOM.FindRedSigma(Sigma,cSigma,wf);
   // for(int i=1; i<nkpoints; i++)
   //   cout<<k[i]<<"\t"<<real(Sigma(i,i))<<"\t"<<imag(Sigma(i,i))<<"\t"<<wf[i]<<"\n";


    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
     p_spectralFunction dom_nonlocal(Mu,n,l,J);
   //  dom_nonlocal.get_ParticleSpectralFunction2();

   vector<double> r(96);
    dom_nonlocal.Sp_rr(Ecm,r);
   Ecm = Ecm+dE;
  
  }
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ////////////////////////////////////////////
   return 0;
}
