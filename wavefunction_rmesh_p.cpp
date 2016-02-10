#include"hartreeFock.h"
#include"boundRspace.h"
// #include"reaction.h"
int main()
{

////////Input Parameters of the potential (fit parameters) /////
std::string parameters_filename="hosfit26Nove.inp";

NuclearParameters Nu = read_nucleus_parameters( "Input/pca40.inp" );

double Ef=-4.7;
int lmax=5;
double z0=20.0;
double zp0;
double A0=40.0;
double tz=0.5;

int type=1;
int mvolume = 4;
int AsyVolume = 1;

double A = 40.0;

if (tz>0) { zp0=1;}
else {zp0=0;}

double ph_gap = 0;
double  rStart = .05;
double  rmax = 12.;
double  ham_pts = 180; 

double  rdelt = rmax / ham_pts;

        // Construct Parameters Object
        Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, zp0 );

        // Construct Potential Object
        pot pottt = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
	pot * pott = &pottt;

        // store coulomb potential in order to compare with




 boundRspace initiate(rmax , ham_pts , Ef, ph_gap , lmax , z0 , zp0 , A0 , pott);
  
 double Elower = -11.61818;
 double Eupper = -9.4;
 double jj = .5;
 int ll = 0;
 int Ifine = 1;
 initiate.searchNonLoc( Elower, Eupper, jj,  ll,  Ifine);
 initiate.exteriorWaveFunct(ll);
 initiate.normalizeWF();


double tol=.05;
double estart=Ef;

///// Making rmesh///

// std::vector<double> rmesh_p= initiate.make_rmesh_point();
 std::vector<double> rmesh= initiate.make_rmesh();


////////////////////////////
//// s and d wave functions////
///////////////////////////


eigen_t waves_s=initiate.find_boundstate(rmesh, estart,1,0,.5,tol);

eigen_t waves_d=initiate.find_boundstate(rmesh, estart,0,2,1.5,tol);

std::ofstream filee("waves/s2andsymp_wave.out");
std::cout<<Elower<<std::endl;
//////////
////////////

double S_s = initiate.sfactor( rmesh, waves_s.first , 0, 0.5, waves_s.second );
double S_d = initiate.sfactor( rmesh, waves_s.first , 2, 1.5, waves_s.second );
// Plottinge Pot and wavefuntionsa and their log//
// and also the derivitaves
//////////////////////o/
double norm_s = 0.;
double norm_d = 0.;
double normAsymp_s = 0.;
double normAsymp_d = 0.;
double Mproton = 931.5;// 937.27;
double normbob= 0.;

double Ed = -waves_d.first;
double Es = -waves_s.first;
double bbetta_d = std::sqrt(2. * Mproton * (40./41.) * Ed /(197. * 197.));
double bbetta_s = std::sqrt(2. * Mproton * (40./41.) * Es /(197. * 197.));
double llambbda = 2.0 * Mproton * (40./41.) * 1.44 * 19. /(197. * 197.);

for(int ii=0;ii<rmesh.size();++ii){
      norm_s += rdelt * std::pow(rmesh[ii] * waves_s.second[ii],2) ;
      norm_d += rdelt * std::pow(rmesh[ii] * waves_d.second[ii],2) ;
      normAsymp_s += rdelt * std::pow((exp(-bbetta_s * rmesh[ii]) * std::pow(rmesh[ii], -llambbda /( 2.0 * bbetta_s))) ,2); 
      normAsymp_d += rdelt * std::pow((exp(-bbetta_d * rmesh[ii]) * std::pow(rmesh[ii], -llambbda /( 2.0 * bbetta_d))) ,2); 
      normbob += rdelt * std::pow(rmesh[ii] * initiate.WaveArray[ii],2);
}
////////
//normalizing wavefunctions at 4.008 to s wave to compare the asymptotic behaviour and energy
//the value of s at 4.008 is .023309//d is .132997 , s_asymp is .000271922 , d_asymp  i s .000336918
/////////
for(int ii=0;ii<rmesh.size();++ii){
      double asymp_s = exp(-bbetta_s * rmesh[ii]) * std::pow(rmesh[ii], -llambbda /( 2.0 * bbetta_s))/(rmesh[ii]); 
      double asymp_d = exp(-bbetta_d * rmesh[ii]) * std::pow(rmesh[ii], -llambbda /( 2.0 * bbetta_d))/(rmesh[ii]); 
//      filee << rmesh[ii] << " " << 2.*std::pow(waves_s.second[ii]/norm_s,2) * S_s
//	    << " "<< 4. * std::pow( waves_d.second[ii]/norm_d ,2) *(.023309/0.132997) * S_d
//            << " " << std::pow( asymp_s/normAsymp_s ,2)  *(.023309/0.000271922) 
//            << " " << std::pow( asymp_d/normAsymp_d ,2)  *(.023309/0.000336918) <<" " <<std::pow( initiate.WaveArray[ii],2) <<std::endl;
      filee << rmesh[ii]
            << " " << 2.0 * std::pow(waves_s.second[ii]/norm_s,2) * S_s 
	    << " "<< 4. * std::pow( waves_d.second[ii]/norm_d ,2) * S_d  
            << " " << 2. * std::pow( asymp_s/normAsymp_s ,2)* (.00158899/(3.10569e-10)) //*(.0038286 /(5.88958e-10)) 
            << " " <<4. * std::pow( asymp_d/normAsymp_d ,2)*(.00542925 /(6.796655e-10)) 
            << " " <<2. * std::pow( initiate.WaveArray[ii]/normbob,2) * (.00158899/.000705115) 
	    <<std::endl;
}
std::cout <<"Es =" << waves_s.first << " " << "Ed= " << waves_d.first <<std::endl;
std::cout << zp0 << " " << tz<<" "<<norm_s<<" "<< norm_d<<" "<< normAsymp_s<<" norm bob=" << normbob<< std::endl;
std::cout <<"S_s =" << S_s << " " << "S_d = " <<S_d  <<std::endl;
std::cout <<"Elower = "<< Elower <<" " << "Eupper = " << Eupper<<std::endl;
filee.close();
return 0;

}
	
