/********************************************************
 *        constants.h
 * DECLARATION OF MATH CONSTANTS
 * AND SOME PHYSICS CONSTANTS
 *
 ********************************************************/
 #ifndef constants_h
 #define constants_h
 #include<cmath> 
 #include<vector>

 using namespace std;
/*********************************************
 *      Numerical CONSTANTS    
 *********************************************/
  double const pi=M_PI;
  double const e_=2.718281828459;
 int  const total_lines = 100;
 double const ln_2=0.6931471805;


/*********************************************
 *        Physical Constants           
 *********************************************/
//  double const M=938.26375;             //Bob's Mass [MeV]
  double const M=938.95;             //Nucleon Mass [MeV]
  double const me=0.511;             //electron mass [MeV]
  double const m_mu=105.6583;        //muon mass     [MeV]
  double const mu_tau=1776.84;       //tau mass         [MeV]
  double const c_=2.99792e+10;      //speed of light    [cm/s]
  double const G_=6.6726e-8;           // Gravitation const [cm3/gs2]
  double const hbarc=197.32858;           //  (hbar*c)       [MeV*fm]
//  double const alpha=0.091701236; //(e^2/(4pi*hbarc)) Structure fine constant
  double const m0=931.5;              //Atomic mass unit [MeV]
  double const kconstant=0.048192;   //  2/hbar^2
  //double const kconstant=0.04822725;   //  2/hbar^2
  double const Mca40 = 39.96259098; // in atomtic units


#endif
