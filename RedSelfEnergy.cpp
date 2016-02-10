/*
                 RedSelfEnergy.cpp


		    ^	^
		      U
		    \___/
    
   -  Reducible Self-Energy definitions

 */
 #include "RedSelfEnergy.hpp"
 using namespace std;
ReducibleSelfEnergy::ReducibleSelfEnergy(int radialQuantumN,
                                         int OrbAngularL,
                                         double TotAngularJ,
                                         double EnergyCM,
					 double ReducedMass,
                                         vector<double> &wvVector,
                                         vector<double> &dwvVector,
                                         vector<double> &FreeProp)
 {
   n = radialQuantumN;
   l = OrbAngularL;
   j = TotAngularJ;
   Ecm = EnergyCM;
   Mu = ReducedMass;
   k  = wvVector;
   dk = dwvVector;
   Go = FreeProp;
 }


/*
   Vdom  must be defined as  MatrixXcd  Vdom(ksize,ksize)

   Outputs:
       R  : DOM  R matrix 
       wf : DOM wave function 
 */
void ReducibleSelfEnergy::FindRMatrix(MatrixXcd &R,
                                      MatrixXcd &cR,
                                      vector<double> &wf )
 { 
   int const ksize=k.size();
   double G0_k2dk;

   MatrixXcd Vdom(ksize,ksize),cVdom(ksize,ksize);
   MatrixXcd U(ksize,ksize),cU(ksize,ksize);
   MatrixXcd U_1(ksize,ksize),cU_1(ksize,ksize);
   MatrixXcd One(ksize,ksize);
   One = MatrixXcd::Identity(ksize,ksize);

   /*  r mesh  */
   double const rmax = 12.;
   int const  rpts = 200;
   
   double rdelt = rmax/rpts;
   vector<double> r, dr;
   dr.assign(rpts,rdelt);
   for(int i =0; i<rpts; ++i)
       r.push_back((i + 0.5)*rdelt );

   /* non-local dom in k-space*/

   dom_nonlocal_r dom_rr(Ecm, n, l, j, r, dr);
   Vdom = dom_rr.dom_k_space( k, wf);

   double norm = 0.;

   for(int ii=0; ii<ksize; ii++) {
        norm += wf[ii] * wf[ii] * k[ii] * k[ii] *dk[ii];
     }
 //  cout<<norm<<endl;

   cVdom =Vdom.conjugate();    // it is evaluated -i eta
 
   for(int ii=0; ii<ksize; ii++)
     {
      wf[ii] /= norm; 
      for(int jj=0; jj<ksize; jj++)
	 {
	   G0_k2dk  = Go[jj]*k[jj]*k[jj]*dk[jj];
	   U(ii,jj) = -G0_k2dk*Vdom(ii,jj);
	  cU(ii,jj) = -G0_k2dk*cVdom(ii,jj);
         }
     }

   /* Finding R-matrix*/
    U  += One;
   cU  += One;

    U_1  =  U.inverse();
   cU_1  = cU.inverse();

    R  =  U_1 * Vdom;
   cR  = cU_1 *cVdom;
 }

/*
   Finding the total reducible Self-Energy
 */
 void ReducibleSelfEnergy::FindRedSigma(MatrixXcd &SelfEnergy,
                                        MatrixXcd &cSelfEnergy,
                                        vector<double> &wf)
 {
   int const ksize=k.size();
   double const ko=k[0];
   complex<double> const irho = complex<double>( 0.0, pi*ko*Mu/hbc2);

   complex<double> denom,cdenom;
   int ii,jj;

   MatrixXcd R(ksize,ksize),cR(ksize,ksize);
   FindRMatrix( R, cR, wf);
  
   denom  = 1.0 + irho*R(0,0);
   cdenom = 1.0 - irho*cR(0,0); 
   for( ii=0; ii<ksize; ii++)
     {
       for( jj=0; jj<ksize; jj++)
	  {
	    SelfEnergy(ii,jj) = R(ii,jj) - irho*R(ii,0)*R(0,jj)/denom;
	    
	    
	    cSelfEnergy(ii,jj) = cR(ii,jj) + 
                                   irho*cR(ii,0)*cR(0,jj)/cdenom;

          }
     }
 }

/*
   Finding the total reducible Self-Energy
 */
 void ReducibleSelfEnergy::FindRedSigma(MatrixXcd &SelfEnergy,
                                        MatrixXcd &cSelfEnergy)
                                        
 {
   int const ksize=k.size();
   double const ko=k[0];
   complex<double> const irho = complex<double>( 0.0, pi*ko*Mu/hbc2);

   complex<double> denom,cdenom;
   int ii,jj;

   MatrixXcd R(ksize,ksize),cR(ksize,ksize);
   FindRMatrix( R, cR);
  
   denom  = 1.0 + irho*R(0,0);
   cdenom = 1.0 - irho*cR(0,0); 
   for( ii=0; ii<ksize; ii++)
     {
       for( jj=0; jj<ksize; jj++)
	  {
	    SelfEnergy(ii,jj) = R(ii,jj) - irho*R(ii,0)*R(0,jj)/denom;
	    
	    
	    cSelfEnergy(ii,jj) = cR(ii,jj) + 
                                   irho*cR(ii,0)*cR(0,jj)/cdenom;

          }
     }
 }

 vector<double> ReducibleSelfEnergy::FindRedSigma(MatrixXcd &SelfEnergy,
                                        MatrixXcd &cSelfEnergy, vector<double> r , vector<double> dr)
                                        
 {
   int const ksize=k.size();
   vector<double> boundwave;
   double const ko=k[0];
   complex<double> const irho = complex<double>( 0.0, pi*ko*Mu/hbc2);

   dom_nonlocal_r dom_rr(Ecm, n, l, j, r, dr);
   complex<double> denom,cdenom;
   int ii,jj;
  //*******************************************/////
  //calculating bound state wf in dom_nonlocal_rspace.cpp
  //*******************************************/////
   boundwave =  dom_rr.dom_r_space_wavefunction(r);

   MatrixXcd R(ksize,ksize),cR(ksize,ksize);
   FindRMatrix( R, cR);
  
   denom  = 1.0 + irho*R(0,0);
   cdenom = 1.0 - irho*cR(0,0); 
   for( ii=0; ii<ksize; ii++)
     {
       for( jj=0; jj<ksize; jj++)
	  {
	    SelfEnergy(ii,jj) = R(ii,jj) - irho*R(ii,0)*R(0,jj)/denom;
	    
	    
	    cSelfEnergy(ii,jj) = cR(ii,jj) + 
                                   irho*cR(ii,0)*cR(0,jj)/cdenom;

          }
     }
 return boundwave;
 }
/**************************************************
       R-matrix function for on-shell calculation
   Then wave function needs not
 **************************************************/

/*
   Vdom  must be defined as  MatrixXcd  Vdom(ksize,ksize)

   Outputs:
       R  : DOM  R matrix 

 */
 void ReducibleSelfEnergy::FindRMatrix(MatrixXcd &R,
                                      MatrixXcd &cR ) 
 { 
   int const ksize=k.size();
   double G0_k2dk;

   MatrixXcd Vdom(ksize,ksize),cVdom(ksize,ksize);
   MatrixXcd U(ksize,ksize),cU(ksize,ksize);
   MatrixXcd U_1(ksize,ksize),cU_1(ksize,ksize);
   MatrixXcd One(ksize,ksize);
   One = MatrixXcd::Identity(ksize,ksize);
   //bound state wf

   /*  r mesh  */
   double const rmaxx =12.;
   int const  rptss = 200;
   double rdeltt = rmaxx/rptss;
   vector<double> rmesht, drmesht;
   drmesht.assign(rptss,rdeltt);
   for(int i =0; i<rptss; ++i)
       rmesht.push_back((i + 0.5)*rdeltt );

   /* non-local dom in k-space*/
   /**
   cout<<"\n";
   cout<<"Before invoquing DOM\n";
   cout<<"Construction values:\n";
   cout<<"Energy ="<<Ecm<<" MeV\n";
   cout<<" Quantum Set "<< n<<"\t"<<l<<"\t"<<j<<"\n";
   cout<<" r-mesh :"<<r[0]<<"\t"<<dr[0]<<"\n";
   cout<<" r-mesh :"<<r[10]<<"\t"<<dr[10]<<"\n";
   cout<<"\n";
   **/
   /***   Testing inputs **/
   dom_nonlocal_r dom_rr(Ecm, n, l, j, rmesht, drmesht);


   Vdom = dom_rr.dom_k_space(k);
   cVdom =Vdom.conjugate();    // it is evaluated -i eta

   for(int ii=0; ii<ksize; ii++)
     {
      for(int jj=0; jj<ksize; jj++)
	 {
	   G0_k2dk  = Go[jj]*k[jj]*k[jj]*dk[jj];
	   U(ii,jj) = -G0_k2dk*Vdom(ii,jj);
	  cU(ii,jj) = -G0_k2dk*cVdom(ii,jj);
         }
     }

    //cout<<"HEYY";
   /* Finding R-matrix*/
    U  += One;
   cU  += One;

    U_1  =  U.inverse();
   cU_1  = cU.inverse();

    R  =  U_1 * Vdom;
   cR  = cU_1 *cVdom;

 }



/*
    Reducible Self-Energy at pole
    If only scattering results are needed on
    shell.
 */
void ReducibleSelfEnergy::SelfEnergyAtPole(complex<double> &S_lj,
					   complex<double> &cS_lj )
{
 int const ksize=k.size();
   int ii,jj;
   double ko=k[0];
   complex<double> Sigma_oo, cSigma_oo;
   complex<double> denom,cdenom;
   complex<double> const irho = complex<double>( 0.0, pi*ko*Mu/hbc2);
   MatrixXcd R(ksize,ksize),cR(ksize,ksize);

   FindRMatrix( R, cR);
   //   S-matrix ...
   Sigma_oo  =  R(0,0)/( 1. + irho*R(0,0) );
   cSigma_oo = cR(0,0)/( 1. - irho*cR(0,0) );

//    string S_file = "diffCrossXn/phaseshiftCompar/E21New/S_matraxNSC.out";
 //   std::ofstream Sfile( S_file.c_str() );
   // for (int kk = 0 ; kk < ksize; ++kk ) {
     //  complex<double> S_matrix = R(kk,kk)/( 1. + irho*R(kk,kk) );
       //Sfile<<k[kk]<<"\t"<< real(S_matrix)<<"\t"<<imag(S_matrix)<<"\n";
   //}
   S_lj  = 1.0 - 2.*irho*Sigma_oo;
   cS_lj = 1.0 + 2.*irho*cSigma_oo; /// Got to be + for unitarity
  
}



