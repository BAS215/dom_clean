/***************************

  Construction of non-local 
  DOM irreducible selfenergy
 -Creating an environment for the
  construction requires several steps.
 -DOM nonlocal in r-space.

 ***************************/

 #include"dom_nonlocal_rspace.hpp"
 using namespace std;
 
  dom_nonlocal_r::dom_nonlocal_r(double Ecm, 
                                int radialQuantumNumber,
                                int OrbitalAngularL, 
                                double TotalAngularJ, 
				vector<double> &r,
                                vector<double> &dr)
  {
   Energy_cm =Ecm; 
   n = radialQuantumNumber; 
   l = OrbitalAngularL; 
   J = TotalAngularJ; 
   rmesh = r;
   drmesh = dr;
  }
 
// MatrixXcd dom_nonlocal_r::dom_r_space(vector<double> &wf)
 
 MatrixXcd dom_nonlocal_r::dom_r_space(cvector_t &Vdomlocal , vector<double> &wf)
 {
   int rsize =rmesh.size(); /// old file it is  rpts
   MatrixXcd Vdom(rsize,rsize);
   
   double rmax =rmesh[rsize-1];
     ////////////////////////////////////////////
     ///////////// DOM environment //////////////
     ////////////////////////////////////////////
//    string input_dir = "Input/";

    int lmax = 5;
  
    ////////
    string parameters_filename = "InputCa40.inp";

    string n_filename = "Input/nca40.inp";

    // Create Nuclear Parameter Objects
    NuclearParameters Nu_n = read_nucleus_parameters( n_filename );
    vector< NuclearParameters > Nu_vec;
    Nu_vec.push_back( Nu_n );

    // Read in DOM parameters
    std::ifstream pfile( parameters_filename.c_str() );
    if ( pfile.is_open() !=1 ) {
        std::cout << "could not open file " << parameters_filename << std::endl;
        std::abort();
    }
    pfile.close();
    pfile.clear();

    // Nucleus Object
    const NuclearParameters &Nu = Nu_vec[0];

   // Construct Parameters Object

    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

    // Construct Potential Object
    pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
    pot *U1 = &U;

      //////////////////////////////////////////
      ////   END OF DOM Environtment Stuff /////
      ///////////////////////////////////////////

    /* CALCULATIONS */
  
    //building V(r,r') and store it in V_r_matrix

 //   U1->setEnergy(Energy_cm);
  //  U1-> setAM(l,J);
    double tolerance =0.05;
 
    boundRspace initiate(rmesh[rsize-1],rsize, Nu.Ef, 0, lmax, Nu.Z, Zp, Nu.A, U1);
    eigen_t waves_function =initiate.find_boundstate(rmesh,Nu.Ef,n,l,J,tolerance);

 // cout <<"ECM "<< Energy_cm<<endl;
  for( int i = 0; i < rsize; ++i ) 
     {
       for( int j = 0; j < rsize; ++j ) 
          {

            U1->setEnergy(Energy_cm);
            U1-> setAM(l,J);
	    Vdom(i,j) =  U1->nonlocalPart( rmesh[i], rmesh[j]) ;
          }
       
          U1->setEnergy(Energy_cm);
          U1-> setAM(l,J);
        //  Vdom(i,i) += U1->localPart(rmesh[i]); 
          
         Vdomlocal[i] = U1->localPart(rmesh[i]); 

         wf[i] = waves_function.second[i];
        // wf[i] = waves_function.second[i]/rmesh[i];
        //cout << rmesh[i]<<" "<<wf[i]<<" "<<Energy_cm<<endl;
     }  
  // Trying the bound state, it might need E<0 ???

   return Vdom;
 }

/***************************************************************
  Output :
    Vdom_kk : DOM nonlocal irreducible self-energy in k-space
    wf_k:  DOM wave function in k-space
 ****************************************************************/
 MatrixXcd dom_nonlocal_r::dom_k_space(vector<double> &k, vector<double> &wf_k)
 {
   int const ksize = k.size();
   int const rsize = rmesh.size();

   MatrixXcd Vdom_kk(ksize,ksize);
   wf_k.resize(ksize);

   /*  DOM Irreducible Self-Energy in r-space*/
   vector<double> wf_r(rsize);

   cvector_t Vdomlocal_r(rsize);
   cvector_t Vdomlocal(rsize);
   MatrixXcd Vdom_rr = dom_r_space(Vdomlocal_r, wf_r);
   

   /*  Transformation to k-space*/
   FourierBesselContainer FB(l,rmesh);
   MatrixXd j_l_kr = FB.all_j_l_kr(k);
   double r2_dr_i, r2_dr_j;
   double wfSum = 0.0;
   complex<double> sumInside = zero;
   complex<double> sumOutside = zero;
   complex<double> sumlocal = zero;
   int ik,jk,ir,jr;

   /*
     Transformation to k-space of the bound wave function is
     performed assuming:
        wf_r =  r*u(r) //From Hossein code.
    */
   //Wave function transformation to k-space  made separately ///

   for(ik=0; ik<ksize; ik++)
      {
       wfSum = 0.0;
       for( ir=0; ir<rsize; ir++)
	  {
            r2_dr_i = rmesh[ir]*rmesh[ir]*drmesh[ir];// Only one r, the other is in wf_r//
           // r2_dr_i = rmesh[ir]*drmesh[ir];// Only one r, the other is in wf_r//
            wfSum +=  r2_dr_i * wf_r[ir] * j_l_kr(ir,ik);
          }
       wf_k[ik] = sqrt(2./pi)*wfSum;
      }

   // Transforming Irreducible Self-Energy
   for(ik=0; ik<ksize; ik++)
     {
       for(jk=0; jk<ksize; jk++)
	  {
	    sumOutside = zero;
            sumlocal = zero;
	    for( ir=0; ir<rsize; ir++)
	       {

               //r2_dr_i = rmesh[ir]*rmesh[ir]*drmesh[ir];
                 r2_dr_i = rmesh[ir]*drmesh[ir];
		 sumInside = zero;

		 for( jr =0; jr<rsize; jr++)
		    {

		    //r2_dr_j = rmesh[jr]*rmesh[jr]*drmesh[jr];
		    r2_dr_j = rmesh[jr]*drmesh[jr];
		     sumInside += r2_dr_j*Vdom_rr(ir,jr)*j_l_kr(jr,jk);
		    // sumInside+= r2_dr_j*Vdom_rr(ir,jr)*bess_mtx(jk,jr);

                    }
                sumlocal +=  rmesh[ir] * rmesh[ir] * drmesh[ir]* j_l_kr(ir,ik) * j_l_kr(ir,jk) * Vdomlocal_r[ir];
		sumOutside +=r2_dr_i*sumInside*j_l_kr(ir,ik);
               }

	    Vdom_kk(ik,jk) = sumOutside + sumlocal;
	
          }
     }


   Vdom_kk = 2.* Vdom_kk /pi;

   return Vdom_kk;
 }



/////////////////// Without calculating bound wave function////////////
/******************************************************************
  Same as previous function dom_r_space, but does not calculate 
  bound wave function.
*******************************************************************/
// MatrixXcd dom_nonlocal_r::dom_r_space()
 MatrixXcd dom_nonlocal_r::dom_r_space( cvector_t &Vdomlocal)
 {
   int rsize =rmesh.size(); /// old file it is  rpts
   MatrixXcd Vdom(rsize,rsize);
   double rmax =rmesh[rsize-1];
     ////////////////////////////////////////////
     ///////////// DOM environment //////////////
     ////////////////////////////////////////////

    string parameters_filename = "Input.inp";

    string n_filename = "Input/nca40.inp";
//    string p_filename = "Input/pca40.inp";


    // Create Nuclear Parameter Objects
    NuclearParameters Nu = read_nucleus_parameters( n_filename );

    // Read in DOM parameters
    std::ifstream pfile( parameters_filename.c_str() );
    if ( pfile.is_open() !=1 ) {
        std::cout << "could not open file " << parameters_filename << std::endl;
        std::abort();
    }
    pfile.close();
    pfile.clear();

    // Nucleus Object

   // Construct Parameters Object
    Parameters p = get_parameters( parameters_filename, Nu.A, Nu.Z, Zp );

    // Construct Potential Object
    pot U = get_bobs_pot2( type, mvolume, AsyVolume, tz, Nu, p );
    pot *U1 = &U;

      //////////////////////////////////////////
      ////   END OF DOM Environtment Stuff /////
      ///////////////////////////////////////////

    /* CALCULATIONS */
  
    //building V(r,r') and store it in V_r_matrix


  for( int i = 0; i < rsize; ++i ) 
     {
       for( int j = 0; j < rsize; ++j ) 
          {
    U1->setEnergy(Energy_cm);
    U1-> setAM(l,J);
	    Vdom(i,j) = U1->nonlocalPart( rmesh[i], rmesh[j]) ;
          }
    U1->setEnergy(Energy_cm);
    U1-> setAM(l,J);
       
      // Vdom(i,i) +=  U1->localPart(rmesh[i]); 
      
       Vdomlocal[i]=U1->localPart(rmesh[i]); 
     }  
  // Trying the bound state, it might need E<0 ???

   return Vdom;
 }


/***************************************************************
  Output :
    Vdom_kk : DOM nonlocal irreducible self-energy in k-space
     No bound wave function is calculated with this function
 ****************************************************************/
 MatrixXcd dom_nonlocal_r::dom_k_space(vector<double> &k)
 {
   int const ksize = k.size();
   int const rsize = rmesh.size();

   MatrixXcd Vdom_kk(ksize,ksize);

   cvector_t Vdomlocal_r(rsize);
 //  cvector_t Vdomlocal(ksize);

   /*  DOM Irreducible Self-Energy in r-space*/

   MatrixXcd Vdom_rr = dom_r_space(Vdomlocal_r );

   /*  Transformation to k-space*/
   FourierBesselContainer FB(l,rmesh);
   MatrixXd j_l_kr = FB.all_j_l_kr(k);
   double r2_dr_i, r2_dr_j;

   complex<double> sumInside = zero;
   complex<double> sumOutside = zero;
   complex<double> sumlocal = zero;
   int ik,jk,ir,jr;

   // Transforming Irreducible Self-Energy
   for(ik=0; ik<ksize; ik++)
     {
       for(jk=0; jk<ksize; jk++)
	  {
	    sumOutside = zero;
            sumlocal = zero;
	    for( ir=0; ir<rsize; ir++)
	       {
                 r2_dr_i = rmesh[ir]*drmesh[ir];
                // r2_dr_i = rmesh[ir]*rmesh[ir]*drmesh[ir];
		 sumInside = zero;
		 for( jr =0; jr<rsize; jr++)
		    {
		      r2_dr_j = rmesh[jr]*drmesh[jr];
		    //  r2_dr_j = rmesh[jr]*rmesh[jr]*drmesh[jr];

		      sumInside+= r2_dr_j*Vdom_rr(ir,jr)*j_l_kr(jr,jk);
                    }

                sumlocal += rmesh[ir]*rmesh[ir]*drmesh[ir]*j_l_kr(ir,ik)*j_l_kr(ir,jk) * Vdomlocal_r[ir];
		sumOutside +=r2_dr_i*sumInside*j_l_kr(ir,ik);
		
               }

	    Vdom_kk(ik,jk) = sumOutside + sumlocal;
	
          }
     }


   Vdom_kk = 2.* Vdom_kk /pi;

   return Vdom_kk;
 }
