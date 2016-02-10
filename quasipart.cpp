//This simply creates the elastic wave functions to input into the ee'p code


#include "quasipart.h"
#include"hartreeFock.h"
#include <complex>

int main(){

   double hbarc = 197.327;

   double Mu = M;   //Mass in MeV
   cout<<Mu<<endl;
   complex <double> M_I(0,1);

   ////////Input Parameters of the potential (fit parameters) /////
   std::string parameters_filename="Input.inp";

   NuclearParameters Nu = read_nucleus_parameters( "Input/pca40.inp" );

   double Ef=Nu.Ef;
   int lmax=20;
   double zp0;
   double tz=0.5;

   int type=1;
   int mvolume = 4;
   int AsyVolume = 1;

   double A = Nu.A;

   double  mu = (A)/((A-1.));

   if (tz>0) { zp0=1;}
   else {zp0=0;}

   double ph_gap = Nu.ph_gap;
   cout<<"ph_gap = "<<Nu.ph_gap<<endl;
   double  rStart = .05;
   double  rmax = 13.51999;
   int  ham_pts = 170; 

   double  rdelt = rmax / (ham_pts-1);
   cout<<"rdelt = "<<rdelt<<endl;

   vector<double> dr;
   dr.assign(ham_pts,rdelt);

   // Construct Parameters Object
   Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, zp0 );
   // Construct Potential Object
   pot pottt = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
   pot * pott = &pottt;
   cout<<"kconstant = "<<pott->kconstant<<endl;

   boundRspace initiate(rmax , ham_pts , Ef, ph_gap , lmax , Nu.Z , zp0 , Nu.A, pott);

   double Elower = -11.61818;
   double Eupper = -9.4;
   double jj = .5;
   int ll = 0;
   int Ifine = 1;
   initiate.searchNonLoc( Elower, Eupper, jj,  ll,  Ifine);
   initiate.exteriorWaveFunct(ll);
   initiate.normalizeWF();


   double tol=.01;
   double estart=Ef;

   ///// Making rmesh///

   std::vector<double> rmesh_p= initiate.make_rmesh_point();
   std::vector<double> rmesh= initiate.make_rmesh();

   eigen_t dom_wf_s = initiate.find_boundstate(rmesh, Ef, 1, 0, 0.5, tol);
   eigen_t dom_wf_d = initiate.find_boundstate(rmesh, Ef, 0, 2, 1.5, tol);

   ofstream dfile("wfdomd32.txt");
   ofstream sfile("wfdoms12.txt");

   double sf = initiate.sfactor(rmesh,dom_wf_d.first,2,1.5,dom_wf_d.second);

   sfile<<rmesh[0]<<" "<<dom_wf_s.second[0]<<endl;
   dfile<<"0.0 0.0"<<endl;
   for(int i=0;i<rmesh.size()-1;i++){
      if(i%2==1){
         dfile<<rmesh[i]<<" "<<dom_wf_d.second[i]<<endl;
         sfile<<rmesh[i]<<" "<<dom_wf_s.second[i]<<endl;
      }
   }

   double rms = 0.0;
   for(int i=0;i<rmesh.size();i++){
      rms += pow(dom_wf_d.second[i]*rmesh[i]*rmesh[i],2) * rdelt;
   }
   rms = sqrt(rms);

   cout<<"rms = "<<rms<<endl;

   cout<<"d3/2 spectroscopic factor = "<<sf<<endl;
   cout<<"d3/2 bound energy = "<<dom_wf_d.first<<endl;
   cout<<"s1/2 bound energy = "<<dom_wf_s.first<<endl;
   cout<<"Fermi Energy = "<<Ef<<endl;

   string title = "pca40";
   string* title0 = &title;

   double E0 = 0.52;

   ofstream rfile("react2.txt");

   double Elab = 100.0;

   if (Elab < 0.) return Elab*A/(1.+A);
   // center of mass velocity in units of c
   double vcm = sqrt(Elab*(Elab+2.*scatterRspace::m0))/(Elab+(1.+A)*
         scatterRspace::m0);
   //gamma factor for this velocity
   double gam = 1./sqrt(1.-std::pow(vcm,2));
   double Ecm = (gam-1.)*(1.+A)*scatterRspace::m0 +
      gam*Elab*(A-1.)*scatterRspace::m0/((A+1.)*
            scatterRspace::m0+Elab);

   //scatterRspace scat(Nu.Z,1,Nu.A,pott,title0);

   //scat.getSmatrix(Ecm,Elab);


   //ofstream eout("elast.txt");
   //ofstream rout("ruth.txt");
   //ifstream din("cross.dat");
   //ofstream dout("cross.txt");
   ////ofstream eout("anal.txt");
   //double ruth,cross,angle;
   //double data,error,theta;
   //int count = 22;

   ////Just needed to run this once to adjust the data by dividing out rutherford
   //for(int i=0;i<count;i++){
      //din >> theta >> data >> error;
      //ruth = scat.Rutherford(theta*3.14159/180);
      //data = data / ruth;
      //error = error / ruth;
      //dout << theta << " " << data << " " << error << endl;
   //}
   //din.close();
   //dout.close();

   //// testing out differential cross section

   //count = 1000;
   //for(int i=1;i<count;i++){
      //angle = i *(3.14159/count/2);
      //cross = scat.DifferentialXsection(angle);
      //ruth = scat.Rutherford(angle);
      //angle = angle * (180.0/3.14159);
      //eout<<angle<<" "<<cross/ruth<<endl;
      //rout<<angle<<" "<<ruth<<endl;
   //}
   //eout.close();

}                                               /* end of main */














