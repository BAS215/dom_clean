/**************************************************
   Spectral Function programme
            --------------

  File: spectralFunction.cpp
  Programmer:Helber Dussan
  Contact:hdussan AT physics.wustl.edu
   
  High level instructions the Reducible 
  Self-Energy calculation using One-body 
  green's function method for nucleon-nucleus 
  scattering. 
  This programme finds the Particle Spectral Function
  Specific for CD-Bonn potential.

   Goal: depletion at positive energies 
---------------------------------------------------
 
   [Energy_cm]= MeV
   [ki]       = fm-1
   [Mu]       = MeV

  High level functions to calculate the 
  Particle Spectral Density

**************************************************/
 #include "spectralFunction.hpp"
 using namespace std;

/***************************************************
  High level function to calculate the 
  Particle Spectral Density
     INPUT PARAMETRES
   Mu : Reduced Mass in MeV
   n  : Princ.Quantum Number
   l   : Orbital angular momentum
   j   : Total angular momentum

  - It does not keep track of each contribution
  - Using spectralFunction()

****************************************************/
void  p_spectralFunction::get_ParticleSpectralFunction()
 {
   double Energy_cm =2.0; // MeV
   double dE = 0.42;       // MeV
   complex<double> strength = zero;
   complex<double> total = zero;
   
   for(int i_e =0; i_e<659; i_e++)
     {
       strength = spectralFunction(Energy_cm);
       Energy_cm+=dE;
       total += strength*dE;
     }
  // cout<<" depletionn =\t"<<total<<"\n";
 }

/***************************************************
  High level function to calculate the 
  Particle Spectral Density
     INPUT PARAMETRES
   Mu : Reduced Mass in MeV
   n  : Princ.Quantum Number
   l   : Orbital angular momentum
   j   : Total angular momentum

  - It keeps track of each contribution
  - Using spectralFunction2()

****************************************************/
void p_spectralFunction::get_ParticleSpectralFunction2()
 {
   double Energy_cm =.2; // MeV
   double dE = .2;       // MeV
   complex<double> strength = zero;
   complex<double> total = zero;
   
   for(int i_e =0; i_e<1000; i_e++)
     {
       vector< complex<double> > p_spectral_f;

       p_spectral_f = spectralFunction2(Energy_cm);
       cout<<Energy_cm<<"\t";
       cout<<real(p_spectral_f[0])<<"\t";
       cout<<imag(p_spectral_f[0])<<"\t";
       cout<<real(p_spectral_f[1])<<"\t";
       cout<<imag(p_spectral_f[1])<<"\t";
       cout<<real(p_spectral_f[2])<<"\t";
       cout<<imag(p_spectral_f[2])<<"\t";
       cout<<real(p_spectral_f[3])<<"\t";
       cout<<imag(p_spectral_f[3])<<"\t";
       cout<<real(p_spectral_f[4])<<"\t";
       cout<<imag(p_spectral_f[4])<<"\t";
       cout<<real(p_spectral_f[5])<<"\t";
       cout<<imag(p_spectral_f[5])<<"\n";


       strength =p_spectral_f[5];
       Energy_cm+=dE;
       total += strength*dE;
     }
   cout<<" depletion =\t"<<total<<"\n";
 }

void p_spectralFunction::savespectral(double Energy_cm) {


 //  complex<double> strength = zero;
 //  complex<double> total = zero;
   
       vector< complex<double> > p_spectral_f;

       p_spectral_f = spectralFunction2(Energy_cm);
       cout<<Energy_cm<<"\t";
       cout<<real(p_spectral_f[0])<<"\t";
       cout<<imag(p_spectral_f[0])<<"\t";
       cout<<real(p_spectral_f[1])<<"\t";
       cout<<imag(p_spectral_f[1])<<"\t";
       cout<<real(p_spectral_f[2])<<"\t";
       cout<<imag(p_spectral_f[2])<<"\t";
       cout<<real(p_spectral_f[3])<<"\t";
       cout<<imag(p_spectral_f[3])<<"\t";
       cout<<real(p_spectral_f[4])<<"\t";
       cout<<imag(p_spectral_f[4])<<"\t";
       cout<<real(p_spectral_f[5])<<"\t";
       cout<<imag(p_spectral_f[5])<<"\n";


      // strength =p_spectral_f[5];
       //Energy_cm+=dE;
      // total += strength*dE;


}

/**************************************************
                   Depletion
 **************************************************/
vector< complex<double> > p_spectralFunction::depletion(vector<double> &k,
							vector<double> &dk,
							vector<double> &Go,
                                                         MatrixXcd &dsigma_plus, 
                                                         MatrixXcd &dsigma_menos,
                                                         vector<double> &wf)
 {
   vector< complex<double> > d_lj;
   double const hbc2 = hbarc*hbarc;
   double const Mu_ko = Mu*k[0]/hbc2;
   complex<double> const i_2pi = complex<double>(0., 1./(2.*pi));
   complex<double> const ipi_Mu_ko_2 = complex<double>(0., pi*Mu_ko*Mu_ko/2.);
   int const ksize = wf.size();

   int ik,jk;
   complex<double> d_0 = zero;
   complex<double> d_1 = zero;
   complex<double> d_2 = zero;
   complex<double> d_3 = zero;
   complex<double> d_4 = zero;
   complex<double> dtemp1 = zero;
   vector< complex<double> > d_int(ksize);
   d_0 = Mu_ko*wf[0]*wf[0];

   double k2_dk_uk_Go;
   for(ik =0 ; ik<ksize; ik++)
     {
      k2_dk_uk_Go=k[ik]*k[ik]*dk[ik]*wf[ik]*Go[ik];
      for( jk=0; jk<ksize; jk++)
	 dtemp1 += k[jk]*k[jk]*dk[jk]*dsigma_menos(ik,jk)*wf[jk]*Go[jk];
      d_1 += k2_dk_uk_Go*dtemp1;
      dtemp1 =zero;

      d_2 += k2_dk_uk_Go*dsigma_plus(0,ik);
      d_3 += k2_dk_uk_Go*dsigma_plus(ik,0);

     }

   d_1 *= i_2pi;
   d_2 *= Mu_ko*wf[0]/2.;
   d_3 *= Mu_ko*wf[0]/2.;
   d_4 = -ipi_Mu_ko_2*wf[0]*wf[0]*dsigma_menos(0,0);

   d_lj.push_back(d_0);
   d_lj.push_back(d_1);
   d_lj.push_back(d_2);
   d_lj.push_back(d_3);
   d_lj.push_back(d_4);

   d_lj.push_back(d_0 + d_1 + d_2 + d_3 + d_4);
   return d_lj;
 }
/**************************************************/
 complex<double> p_spectralFunction::spectralFunction(double Energy_cm)
 {
   complex<double> depletion;
    /*************************
      Construction of k mesh
     *************************/
    double ko;
    int nkpoints;
    vector<double> ki,wi,Go_i;
    kgrid k_points(Energy_cm,Mu);
    k_points.makeKmesh(ki,wi);  
    nkpoints = ki.size();
    k_points.getPropagator(ki,Go_i);
    ko=k_points.getPole(); 

    /*
       Reducible Self-Energy
     */
    
    ReducibleSelfEnergy cdbSelfE(n,l,j,Energy_cm,Mu,ki,wi,Go_i);
    MatrixXcd sigmakk(nkpoints,nkpoints);
    MatrixXcd csigmakk(nkpoints,nkpoints);
    vector<double> wf_k;
    cdbSelfE.FindRedSigma(sigmakk,csigmakk,wf_k);
    
    /******************************
      Overlap function projected on 
         bound state function
     *******************************/
    /*    
      Find depletion using bound wave function here!!!!
    
    */
    cout<<"\t"<<Energy_cm<<"\t";
    cout<<real(depletion)<<"\t";  
    cout<<imag(depletion)<<"\n";

    return depletion;     
 }

/************************************************************

    Keeping track of each contribution

 ************************************************************/
 vector< complex<double> > p_spectralFunction::spectralFunction2(double Energy_cm) 
 {
   vector< complex<double> > s_p;

    /*************************
      Construction of k mesh
     *************************/
    double ko;
    int nkpoints;
    vector<double> ki,wi,Go_i;
    kgrid k_points(Energy_cm,Mu);
    k_points.makeKmesh(ki,wi);  
    nkpoints = ki.size();
    k_points.getPropagator(ki,Go_i);
    ko=k_points.getPole(); 

    /*
       Reducible Self-Energy
     */
    
    ReducibleSelfEnergy domSelfE(n,l,j,Energy_cm,Mu,ki,wi,Go_i);
    MatrixXcd sigmakk(nkpoints,nkpoints);
    MatrixXcd csigmakk(nkpoints,nkpoints);
    vector<double> wf_k;
    domSelfE.FindRedSigma(sigmakk,csigmakk,wf_k);

/*
    for (int i = 0;i<nkpoints;++i){
         for (int ii = 0;ii<nkpoints;++ii){
            // cout<<ki[i]<<" " <<ki[ii]<<" "<<real(sigmakk(i,ii))<<endl; 
         //    cout<<real(csigmakk(i,ii))<<" "; 
             } 
        // cout<<endl;
    }
*/

    /******************************
       Overlap function projected on 
         bound state function
    *******************************/
    MatrixXcd splus = sigmakk+ csigmakk; 
    MatrixXcd sminus = sigmakk - csigmakk; 
//    MatrixXcd splus =  sigmakk; 
 //   MatrixXcd sminus = sigmakk; 
    s_p = depletion(ki,wi,Go_i,splus,sminus,wf_k);

    return s_p;
 }

/******************************************************
    Free Particle Spectral Function in position space 
 ******************************************************/
 MatrixXcd p_spectralFunction::Sp0_rr(double ko,MatrixXd &J_kr)
 {
   int const rsize = J_kr.cols();
   double const rho = Mu*ko/hbc2;
   int i,j;
   double j_l_kr_j_l_kr;

   MatrixXcd Sp0(rsize,rsize);

   for(i=0; i<rsize; i++)
     {
       for(j=0; j<rsize; j++)
	   Sp0(i,j) = J_kr(0,i)*J_kr(0,j);
         
     }

   Sp0 *= 2.*rho/pi;
   return Sp0;
 }
/******************************************************
    Contrib. 1 to the Particle Spectral Function 
    in position space (1st correlated)
 ******************************************************/
  MatrixXcd p_spectralFunction::Sp1_rr( int rsize,
                                        vector<double> &dk,
                     	      	        vector<double> &k,
				        vector<double> &G0,
				        MatrixXcd &sigmakkMenos,
				        MatrixXd &J_kr)
 {
   int const ksize =k.size();
   complex<double> const i_pi2 = complex<double> (0.,1/(pi*pi));
   MatrixXcd S1_rr(rsize,rsize);
   complex<double> sumInside = zero;
   complex<double> sumOutside = zero;
   vector<double> k2_dk_G0(ksize);
   int ir,jr,ik,jk;
   for(ik=0; ik<ksize;ik++) k2_dk_G0[ik]=k[ik]*k[ik]*dk[ik]*G0[ik];

   for(ir=0; ir<rsize;ir++)
      {
       for(jr=0; jr<rsize;jr++)
          {
	   sumOutside = zero;
           for(ik=0; ik<ksize;ik++)   
	      {       
	       sumInside = zero;
               for(jk=0; jk<ksize;jk++)
	          {
		   sumInside += k2_dk_G0[jk]*sigmakkMenos(ik,jk)*J_kr(jk,jr);          
                  }
	       sumOutside+= k2_dk_G0[ik]*sumInside*J_kr(ik,ir);
	      }

	    S1_rr(ir,jr) = sumOutside;
           }
     }
   S1_rr *= i_pi2;
   return S1_rr;
 }

/******************************************************
    Contrib. 2 to the Particle Spectral Function 
    in position space (2st correlated)
 ******************************************************/
  MatrixXcd p_spectralFunction::Sp2_rr( int rsize,
                                        vector<double> &dk,
                     	      	        vector<double> &k,
				        vector<double> &G0,
				        MatrixXcd &sigmakkPlus,
				        MatrixXd &J_kr)
  {
    MatrixXcd S2_rr(rsize,rsize);
    double const rho = Mu*k[0]/hbc2;
    int const ksize = k.size();
    int ir,jr,jk;
    complex<double> sumInside = zero;
    double k2_dk_G0;
    for( ir=0; ir<rsize; ir++)
       { 
	 for( jr=0; jr<rsize; jr++)
	    {
              sumInside = zero;
	      for( jk = 0; jk<ksize; jk++ )
		 {
		   k2_dk_G0 =k[jk]*k[jk]*dk[jk]*G0[jk];
		   sumInside += k2_dk_G0*sigmakkPlus(0,jk)*J_kr(jk,jr);
                 }
	      S2_rr(ir,jr) = J_kr(0,ir)*sumInside;
            }
       }
  
    S2_rr  *= rho/pi;
 return S2_rr;
}

/******************************************************
    Contrib. 3 to the Particle Spectral Function 
    in position space (3rd correlated)
 ******************************************************/
  MatrixXcd p_spectralFunction::Sp3_rr( int rsize,
                                        vector<double> &dk,
                     	      	        vector<double> &k,
				        vector<double> &G0,
				        MatrixXcd &sigmakkPlus,
				        MatrixXd &J_kr)
  {   
    MatrixXcd S3_rr(rsize,rsize);
    double const rho = Mu*k[0]/hbc2;
    int const ksize = k.size();
    int ir,jr,ik;
    complex<double> sumInside = zero;
    double k2_dk_G0;
    for( ir=0; ir<rsize; ir++)
       { 
	 for( jr=0; jr<rsize; jr++)
	    {
              sumInside = zero;
	      for( ik = 0; ik<ksize; ik++ )
		 {
		   k2_dk_G0 =k[ik]*k[ik]*dk[ik]*G0[ik];
		   sumInside += k2_dk_G0*sigmakkPlus(ik,0)*J_kr(ik,ir);
                 }
	      S3_rr(ir,jr) = J_kr(0,jr)*sumInside;
            }
       }

    S3_rr  *= rho/pi;
    return S3_rr;

  }

/******************************************************
    Contrib. 4 to the Particle Spectral Function 
    in position space (4th correlated)
 ******************************************************/
 MatrixXcd p_spectralFunction::Sp4_rr(double ko,
				      MatrixXcd &sigmakkMenos,
                                      MatrixXd &J_kr)
 {
   int const rsize = J_kr.cols();
   double const rho = Mu*ko/hbc2;
   complex<double> const i_ = complex<double> (0.0,1.0);
   int i,j;
   double j_l_kr_j_l_kr;
   MatrixXcd Sp4(rsize,rsize);
   for(i=0; i<rsize; i++)
      {
       for(j=0; j<rsize; j++)
	   Sp4(i,j) = J_kr(0,i)*sigmakkMenos(0,0)*J_kr(0,j);        
      }

   Sp4 *= i_*rho*rho;
   return Sp4;
 }

/*************************************************************************
 * the 2 contributions in the wave function for E>0
 *
 *
 */
   vector< complex<double> > p_spectralFunction::Sp0wave_r( int rsize,
                                        vector<double> &dk,
                     	      	        vector<double> &k,
				        vector<double> &G0,
				        MatrixXcd &sigmakkPlus,
				        MatrixXd &J_kr)
  {
    int const ksize = k.size();
    int ir,jr,jk;
    double k2_dk_G0;
    double const rho = Mu*k[0]/hbc2;
   vector< complex <double> > Sp0;
    for( ir=0; ir<rsize; ir++)
       { 
        complex<double> sumInside = zero;
        for( jk = 0; jk<ksize; jk++ )
	   {
	   k2_dk_G0 =k[jk]*k[jk]*dk[jk]*G0[jk];
	   sumInside += k2_dk_G0*sigmakkPlus(0,jk)*J_kr(jk,ir);
        }
    Sp0.push_back(sumInside);
   }
    return Sp0;
  }

   vector< complex<double> > p_spectralFunction::Sp1wave_r( int rsize,
                                        vector<double> &dk,
                     	      	        vector<double> &k,
				        vector<double> &G0,
				        MatrixXcd &sigmakkPlus,
				        MatrixXd &J_kr)
  {
    
    int const ksize = k.size();
    int ir,jr,jk;
    complex<double> sumInside = zero;
    double k2_dk_G0;
    double const rho = Mu*k[0]/hbc2;
    complex<double> coefff = -complex<double> (0. , 1.)*pi*rho;
   vector< complex<double> > Sp1;
        for( ir = 0; ir<rsize; ir++ )
	   {
	    sumInside = (coefff*sigmakkPlus(0,0)*J_kr(0,ir));
            cout<<"sigmazero   " <<sigmakkPlus(0,0)<<endl;
	    Sp1.push_back(sumInside);

        }
    return Sp1;
  }
//****************************************************************
//
//
//
/******************************************************
 *                                                    *
 *   Full Correlated Particle Spectral Function       *
 *           in position space                        *
 *                                                    *
 *   Input:                                           *
 *       Energy_cm := Energy in the cm frame [MeV]    *
 *                                                    *
 *   Output:                                          *
 *        r      := radial position [fm]              *
 *        Sp_nlj := Particle Spectral Function        *
 *                  [MeV-1 fm-3]                      *
 *                                                    *
 ******************************************************/
// MatrixXcd p_spectralFunction::Sp_rr(double Energy_cm,vector<double> &r)
 MatrixXcd p_spectralFunction::Sp_rr(double Energy_cm,vector<double> r,vector<double> dr)
 {
   MatrixXcd Sp_nlj;
   /*************************
      Construction of k mesh
    *************************/
    double ko;
    int nkpoints;
    vector<double> k,dk,Go;
    kgrid k_points(Energy_cm,Mu);
    k_points.makeKmesh(k,dk);  
    nkpoints = k.size();
    k_points.getPropagator(k,Go);
    ko=k_points.getPole(); 

    /*
       Reducible Self-Energy
     */
    
    ReducibleSelfEnergy cdbSelfE(n,l,j,Energy_cm,Mu,k,dk,Go);
    MatrixXcd sigmakk(nkpoints,nkpoints);
    MatrixXcd csigmakk(nkpoints,nkpoints);
    vector<double> wf_k;    
   // cdbSelfE.FindRedSigma(sigmakk,csigmakk,wf_k);
    vector<double> wavebound = cdbSelfE.FindRedSigma(sigmakk,csigmakk,r,dr);


    /*******************************
     Fourier-Bessel Transformation
      to position space
         r from 0 to 12 fm
     *******************************/

    MatrixXcd sigmakkPlus = sigmakk + csigmakk;
    MatrixXcd sigmakkMenos = sigmakk- csigmakk;
  //  MatrixXcd sigmakkPlus =  sigmakk;
  //  MatrixXcd sigmakkMenos = sigmakk;
    int rsize = r.size();
   // vector<double> dr(rsize);
    
   //for( int i = 0; i < 200; ++i ) {

   //     dr.push_back( 10./200. );
   //     r.push_back( ( i + 0.5 ) * dr );
   // }
   // GausLeg(0.,12.,r,dr);    

    FourierBesselContainer to_rspace(l,k);
    MatrixXd J_kr = to_rspace.all_j_l_kr(r);
    
    /*  Keeping track of each contribution */ 
    MatrixXcd Sp0,Sp1,Sp2,Sp3,Sp4;
    vector < complex <double> > Sp0w , Sp1w; 

   // Sp0w = Sp0wave_r(rsize,dk,k,Go,sigmakkPlus,J_kr);
   // Sp1w = Sp1wave_r(rsize,dk,k,Go,sigmakkPlus,J_kr);

    Sp0 = Sp0_rr(ko,J_kr);
    Sp1 = Sp1_rr(rsize,dk,k,Go,sigmakkMenos,J_kr);
    Sp2 = Sp2_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
    Sp3 = Sp3_rr(rsize,dk,k,Go,sigmakkPlus,J_kr);
    Sp4 = Sp4_rr(ko,sigmakkMenos,J_kr);

    Sp_nlj = Sp0 + Sp1 + Sp2 + Sp3 - Sp4;

  //  double eta = 0.5353;
  //  double onePluseta_2 = (1+eta)*(1+eta);      
  //  double cos2delta = -.230;
  //  double sin2delta = 0.694;
  //  double eta = .1066;// .5380;//.5007;//.4219;//.4432;//.5040;//0.5381;
  //  double onePluseta_2 = (1+eta)*(1+eta);      
  //  double cos2delta = .5899;//-.5194;//-.779;//.0672;//.8182;// 0.879;// -.5195;//-.230;
  //  double sin2delta = .8075;//+.8544;//.6277; -.998;// -.5749;//0.476;// .8545;//0.694;
    double sumIn;
    double sumOut;    
    sumOut = 0;
    double norm=0;
    for(int ir=0; ir<rsize; ir++)
       {
//	 cout<<Energy_cm<<"\t";
/*	 cout<<r[ir]<<"\t";
         cout<<r[ir]*r[ir]*real(Sp0(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp0(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*real(Sp1(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp1(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*real(Sp2(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp2(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*real(Sp3(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp3(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*real(Sp4(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp4(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*real(Sp_nlj(ir,ir))<<" ";
	 cout<<r[ir]*r[ir]*imag(Sp_nlj(ir,ir))<<"\n ";
*/
//       cout<<pi*hbc2*r[ir]*r[ir]*real(Sp_nlj(ir,ir))/(2.*ko*Mu)<<" ";
	 
//       double Snlj = pi*hbc2*r[ir]*r[ir]*real(Sp_nlj(ir,ir))/(2.*Mu);
//	 cout<<Snlj<<"\n";
//	 cout<<r[ir]*real(Sp0w[ir])<<"\t";
//	 cout<<r[ir]*imag(Sp0w[ir])<<"\t";
//	 cout<<r[ir]*real(Sp1w[ir])<<"\t";
//	 cout<<r[ir]*imag(Sp1w[ir])<<"\n";
//	 cout<<  r[ir] * J_kr(0,ir) * r[ir] * J_kr(0,ir) <<"\t";
// 	 cout<< .25* (onePluseta_2 - 2.0 * eta *(1. + (cos(2.*ko*r[ir])*cos2delta-sin2delta*sin(2.*ko*r[ir]))))/(ko * ko) <<"\t";
         double wave_real = (J_kr(0,ir)+real(Sp0w[ir])+real(Sp1w[ir]));
         double wave_imag = (imag(Sp0w[ir])+imag(Sp1w[ir]));
//         double r2_wave_2 = ko * r[ir] * r[ir] * (wave_real * wave_real+ wave_imag * wave_imag);
//	 cout<< r2_wave_2 <<" "<<endl;
//	 cout<< Snlj <<" ";
	 norm += dr[ir]*r[ir]*r[ir]*wavebound[ir]* wavebound[ir];
         sumIn=0;
	 for (int jr=0; jr<rsize; jr++)
             {

	       sumIn = sumIn + dr[jr]*r[jr]*r[jr]*real(Sp_nlj(ir,jr))*wavebound[jr];
            //   cout<<"Sum In = " <<sumIn<<endl;
     	     }
         sumOut +=  dr[ir] * r[ir] *r[ir] * wavebound[ir] * sumIn; 
//         cout<<r[ir]<<" "<<wavebound[ir]<<endl;
 //        cout<<sumOut<<endl;
       }
       
//	 cout<<"\n";
   cout<<Energy_cm<<"\t"<<sumOut<<"\t"<<sumOut*pi*hbc2/(2.0*Mu)<<"\t"<<sumOut*pi*hbc2/(2.0*Mu*ko)<<" "<<norm<<endl;
   return Sp_nlj;
 }
