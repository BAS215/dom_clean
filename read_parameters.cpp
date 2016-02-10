
#include "read_parameters.h"

Raw_Parameters read_par_from_file( const std::string &parname ) {

    std::ifstream file(parname.c_str());

    // if one cannot open file quit
    if (!file.is_open()) 
    {
        std::cout << "couldn't open data file " << parname << std::endl;
        std::abort();
    }

    int mvolume;
    int iFermi;
    file >> mvolume;
    file >> iFermi;

    static int const TotPara=101;
    int varied[TotPara];
    double scaling[TotPara];
    double allPara[TotPara];
    int squared[TotPara];  

    std::string variable;
    for (int i=0;i<TotPara;i++) {

        file >> allPara[i] >> varied[i] >> squared[i] >> scaling[i] >> variable;
        if (squared[i]) allPara[i] = std::sqrt(allPara[i]); 

        //std::cout << variable << ": " << allPara[i] << std::endl;
    }
    file.close();
    file.clear();

    Raw_Parameters p;

    int index = 0;

    p.rc0 = allPara[index++];
    p.rc1 = allPara[index++];
    p.VHFvol = allPara[index++];
    p.VHFsur = allPara[index++];
    p.beta_nl_R0 = allPara[index++];
    p.beta_nl_R1 = allPara[index++];
    p.beta_nl_I0 = allPara[index++];
    p.beta_nl_I1 = allPara[index++];
    p.beta_nl_I0_sur = allPara[index++];
    p.beta_nl_I1_sur = allPara[index++];
    p.AHF = allPara[index++];
    p.rHF0 = allPara[index++];
    p.rHF0s = allPara[index++];
    p.rHF1 = allPara[index++];
    p.aHF0 = allPara[index++];
    p.aHF0s = allPara[index++];
    p.aHF1 = allPara[index++];
    p.fGap_A = allPara[index++];
    p.fGap_B = allPara[index++];
    p.BsurfaceA = allPara[index++];
    p.CsurfaceA = allPara[index++];
    p.DsurfaceA = allPara[index++];
    p.Bsurface = allPara[index++];
    p.Csurface = allPara[index++];
    p.Dsurface = allPara[index++];

    p.Asurface9P = allPara[index++];
    p.Asurface9N = allPara[index++];
    p.Asurface40P_A = allPara[index++];
    p.Asurface40N_A = allPara[index++];
    p.Asurface40P_B = allPara[index++];
    p.Asurface40N_B = allPara[index++];
    p.Asurface42P = allPara[index++];
    p.Asurface44P = allPara[index++];
    p.Asurface48P = allPara[index++];
    p.Asurface48N = allPara[index++];
    p.Asurface50P = allPara[index++];
    p.Asurface52P = allPara[index++];
    p.Asurface52N = allPara[index++];
    p.Asurface54P = allPara[index++];
    p.Asurface54N = allPara[index++];
    p.Asurface58P = allPara[index++];
    p.Asurface58N = allPara[index++];
    p.Asurface60P = allPara[index++];
    p.Asurface60N = allPara[index++];
    p.Asurface62P = allPara[index++];
    p.Asurface62N = allPara[index++];
    p.Asurface64P = allPara[index++];
    p.Asurface64N = allPara[index++];
    p.Asurface88P = allPara[index++];
    p.Asurface90P = allPara[index++];
    p.Asurface92P = allPara[index++];
    p.Asurface92N = allPara[index++];
    p.Asurface112P = allPara[index++];
    p.Asurface114P = allPara[index++];
    p.Asurface116P = allPara[index++];
    p.Asurface116N = allPara[index++];
    p.Asurface118P = allPara[index++];
    p.Asurface118N = allPara[index++];
    p.Asurface120P = allPara[index++];
    p.Asurface120N = allPara[index++];
    p.Asurface122P = allPara[index++];
    p.Asurface124P = allPara[index++];
    p.Asurface124N = allPara[index++];
    p.Asurface206P = allPara[index++];
    p.Asurface208P = allPara[index++];
    p.Asurface208N = allPara[index++];

    p.rsurface0Above = allPara[index++];
    p.rsurface0Below = allPara[index++];

    p.rsurface1 = allPara[index++];
    p.asurfaceAbove = allPara[index++];
    p.asurfaceBelow = allPara[index++];
//NEW
    p.EpvolumeAbove = pow(allPara[index++],2);
    p.EpvolumeBelow = pow(allPara[index++],2);
//
    p.Avolume0_A = allPara[index++];
    p.Avolume0_B = allPara[index++];
//
    p.Avolume1 = allPara[index++];
//NEW
    p.Bvolume0Above = allPara[index++];
    p.Bvolume0Below = allPara[index++];
//
    p.Bvolume1 = allPara[index++];
//NEW
    p.rvolume0Above = allPara[index++];
    p.rvolume0Below = allPara[index++];
//
    p.rvolume1 = allPara[index++];
    p.deltaRvolume = allPara[index++];
    p.expRvolume = allPara[index++];
//NEW
    p.avolume0Above = allPara[index++];
    p.avolume0Below = allPara[index++];
//
    p.avolume1 = allPara[index++];
    p.Ea_above = allPara[index++];
    p.Ea_below = allPara[index++];
    p.alphaOverA = allPara[index++];
    p.Vso = allPara[index++];
    p.VsoNZ = allPara[index++];
    p.rso0 = allPara[index++];
    p.rso1 = allPara[index++];
    p.aso0 = allPara[index++];
    p.aso1 = allPara[index++];
    p.AWso = allPara[index++];
    p.BWso = allPara[index++];
    p.V_wine = allPara[index++];
    p.R_wine = allPara[index++]; 
    p.rho_wine = allPara[index++]; 
   return p;
}

/*
std::vector<double> 
Parameters_to_vector( const Raw_Parameters &rp ) {

    std::vector<double> vec;

    vec.push_back( p.rc0 );
    vec.push_back( p.rc1 );
    vec.push_back( p.VHFvol );
    vec.push_back( p.VHFsur );
    vec.push_back( 
    
}
*/

Parameters
get_parameters( const std::string &filename, double A, double Z, double Zp ) {

    //asymmetry = (N-Z)/A
    double asymmetry;
    if( Zp == 1 ) asymmetry = ( A - 2 * Z ) / A;
    else asymmetry = - ( A - 2 * Z ) / A;

    // raw parameters
    Raw_Parameters p = read_par_from_file( filename );
    
    // modified parameters
    Parameters mp;

    mp.Rc = std::pow( A, 1./3. ) * p.rc0 + p.rc1;
    mp.VHFvol = p.VHFvol;
    mp.VHFsur = p.VHFsur;
    mp.beta_nl_R0 = p.beta_nl_R0;
    mp.beta_nl_R1 = p.beta_nl_R1;
    mp.beta_nl_I0 = p.beta_nl_I0;
    mp.beta_nl_I1 = p.beta_nl_I1;
    mp.beta_nl_I0_sur = p.beta_nl_I0_sur;
    mp.beta_nl_I1_sur = p.beta_nl_I1_sur;
    mp.AHF = p.AHF;
    mp.RHF = p.rHF0 * std::pow(A,1./3.) + p.rHF1;
    mp.aHF = p.aHF0 + p.aHF1 / std::pow(A,1./3.);
    mp.RHFs = p.rHF0s * std::pow(A,1./3.) +p.rHF1;
    mp.aHFs = p.aHF0s + p.aHF1 / std::pow(A,1./3.);
    mp.fGap_A = p.fGap_A;
    mp.fGap_B = p.fGap_B;
    mp.BsurfaceA = p.BsurfaceA;
    mp.CsurfaceA = p.CsurfaceA;
    mp.DsurfaceA = p.DsurfaceA;
    mp.Bsurface = p.Bsurface;
    mp.Csurface = p.Csurface;
    mp.Dsurface = p.Dsurface;
    
    if (Zp == 1 && A == 9) {

        mp.AsurfaceAbove = p.Asurface9P;
        mp.AsurfaceBelow = p.Asurface9P;
    }
    else if (Zp == 0 && A == 9) {

        mp.AsurfaceAbove = p.Asurface9N;
        mp.AsurfaceBelow = p.Asurface9N;
    }
    else if (Zp == 1 && A == 40) {
    
        mp.AsurfaceAbove = p.Asurface40P_A;
        mp.AsurfaceBelow = p.Asurface40P_B;
    }
    else if (Zp == 0 && A == 40) {

        mp.AsurfaceAbove = p.Asurface40N_A;
        mp.AsurfaceBelow = p.Asurface40N_B;
    }
    else if (Zp == 1 && A == 48) {

        mp.AsurfaceAbove = p.Asurface48P;
        mp.AsurfaceBelow = p.Asurface48P;
    }
    else if (Zp == 0 && A == 48) {

        mp.AsurfaceAbove = p.Asurface48N;
        mp.AsurfaceBelow = p.Asurface48N;
    }
    else if (Zp == 1 && A == 42) {

        mp.AsurfaceAbove = p.Asurface42P;
        mp.AsurfaceBelow = p.Asurface42P;
    }
    else if (Zp == 1 && A == 44) {

        mp.AsurfaceAbove = p.Asurface44P;
        mp.AsurfaceBelow = p.Asurface44P;
    }
    else if (Zp == 1 && A == 50) {

        mp.AsurfaceAbove = p.Asurface50P;
        mp.AsurfaceBelow = p.Asurface50P;
    }
    else if (Zp == 1 && A == 52) {

        mp.AsurfaceAbove = p.Asurface52P;
        mp.AsurfaceBelow = p.Asurface52P;
    }
    else if (Zp == 0 && A == 52) {

        mp.AsurfaceAbove = p.Asurface52N;
        mp.AsurfaceBelow = p.Asurface52N;
    }
    else if (Zp == 1 && A == 58) {

        mp.AsurfaceAbove = p.Asurface58P;
        mp.AsurfaceBelow = p.Asurface58P;
    }
    else if (Zp == 1 && A == 60) {

        mp.AsurfaceAbove = p.Asurface60P;
        mp.AsurfaceBelow = p.Asurface60P;
    }
    else if (Zp == 1 && A == 62) {

        mp.AsurfaceAbove = p.Asurface62P;
        mp.AsurfaceBelow = p.Asurface62P;
    }
    else if (Zp == 1 && A == 64) {

        mp.AsurfaceAbove = p.Asurface64P;
        mp.AsurfaceBelow = p.Asurface64P;
    }
    else if (Zp == 0 && A == 58) {

        mp.AsurfaceAbove = p.Asurface58N;
        mp.AsurfaceBelow = p.Asurface58N;
    }
    else if (Zp == 0 && A == 60) {

        mp.AsurfaceAbove = p.Asurface60N;
        mp.AsurfaceBelow = p.Asurface60N;
    }
    else if (Zp == 0 && A == 62) {

        mp.AsurfaceAbove = p.Asurface62N;
        mp.AsurfaceBelow = p.Asurface62N;
    }
    else if (Zp == 0 && A == 64) {

        mp.AsurfaceAbove = p.Asurface64N;
        mp.AsurfaceBelow = p.Asurface64N;
    }
    else if (Zp == 1 && A == 54) {

        mp.AsurfaceAbove = p.Asurface54P;
        mp.AsurfaceBelow = p.Asurface54P;
    }
    else if (Zp == 0 && A == 54) {

        mp.AsurfaceAbove = p.Asurface54N;
        mp.AsurfaceBelow = p.Asurface54N;
    }
    else if (Zp == 1 && A == 88) {

        mp.AsurfaceAbove = p.Asurface88P;
        mp.AsurfaceBelow = p.Asurface88P;
    }
    else if (Zp == 1 && A == 90) {

        mp.AsurfaceAbove = p.Asurface90P;
        mp.AsurfaceBelow = p.Asurface90P;
    }
    else if (Zp == 1 && A == 92) {

        mp.AsurfaceAbove = p.Asurface92P;
        mp.AsurfaceBelow = p.Asurface92P;
    }
    else if (Zp == 0 && A == 92) {

        mp.AsurfaceAbove = p.Asurface92N;
        mp.AsurfaceBelow = p.Asurface92N;
    }
    else if (Zp == 1 && A == 112) {

        mp.AsurfaceAbove = p.Asurface112P;
        mp.AsurfaceBelow = p.Asurface112P;
    }
    else if (Zp == 1 && A == 114) {

        mp.AsurfaceAbove = p.Asurface114P;
        mp.AsurfaceBelow = p.Asurface114P;
    }
    else if (Zp == 1 && A == 116) {

        mp.AsurfaceAbove = p.Asurface116P;
        mp.AsurfaceBelow = p.Asurface116P;
    }
    else if (Zp == 0 && A == 116) {

        mp.AsurfaceAbove = p.Asurface116N;
        mp.AsurfaceBelow = p.Asurface116N;
    }
    else if (Zp == 1 && A == 118) {

        mp.AsurfaceAbove = p.Asurface118P;
        mp.AsurfaceBelow = p.Asurface118P;
    }
    else if (Zp == 0 && A == 118) {
        
        mp.AsurfaceAbove = p.Asurface118N;
        mp.AsurfaceBelow = p.Asurface118N;
    }
    else if (Zp == 1 && A == 120) {

        mp.AsurfaceAbove = p.Asurface120P;
        mp.AsurfaceBelow = p.Asurface120P;
    }
    else if (Zp == 0 && A == 120) {

        mp.AsurfaceAbove = p.Asurface120N;
        mp.AsurfaceBelow = p.Asurface120N;
    }
    else if (Zp == 1 && A == 122) {

        mp.AsurfaceAbove = p.Asurface122P;
        mp.AsurfaceBelow = p.Asurface122P;
    }
    else if (Zp == 1 && A == 124) {

        mp.AsurfaceAbove = p.Asurface124P;
        mp.AsurfaceBelow = p.Asurface124P;
    }
    else if (Zp == 0 && A == 124) {

        mp.AsurfaceAbove = p.Asurface124N;
        mp.AsurfaceBelow = p.Asurface124N;
    }
    else if ( Zp == 1 && A == 206) {

        mp.AsurfaceAbove = p.Asurface206P;
        mp.AsurfaceBelow = p.Asurface206P;
    }
    else if ( Zp == 1 && A == 208) {

        mp.AsurfaceAbove = p.Asurface208P;
        mp.AsurfaceBelow = p.Asurface208P;
    }
    else if ( Zp == 0 && A == 208) {

        mp.AsurfaceAbove = p.Asurface208N;
        mp.AsurfaceBelow = p.Asurface208N;
    }
    else {
        std::cout << "reaction not implemented " << std::endl;
        std::cout << Zp << " " << A << std::endl;
	    abort();
    }

     mp.RsurfaceAbove = std::pow(A,1./3.) * p.rsurface0Above + p.rsurface1;
     mp.RsurfaceBelow = std::pow(A,1./3.) * p.rsurface0Below + p.rsurface1;

     mp.EpvolumeAbove = p.EpvolumeAbove;
     mp.EpvolumeBelow = p.EpvolumeBelow;

     mp.AvolumeAbove = p.Avolume0_A + asymmetry * p.Avolume1;
     mp.AvolumeBelow = p.Avolume0_B + asymmetry * p.Avolume1;

     mp.BvolumeAbove = p.Bvolume0Above + asymmetry * p.Bvolume1;
     mp.BvolumeBelow = p.Bvolume0Below + asymmetry * p.Bvolume1;

     mp.RvolumeAbove = p.rvolume0Above * std::pow(A,1./3.) + p.rvolume1;
     mp.RvolumeBelow = p.rvolume0Below * std::pow(A,1./3.) + p.rvolume1;

     mp.deltaRvolume = p.deltaRvolume; 
     mp.expRvolume = std::pow( p.expRvolume, 2 ); //FIXME why squared?

     mp.avolumeAbove = p.avolume0Above + p.avolume1 / std::pow(A,1./3.);
     mp.avolumeBelow = p.avolume0Below + p.avolume1 / std::pow(A,1./3.);

    // mp.asurfaceAbove = mp.avolumeAbove;
   //  mp.asurfaceBelow = mp.avolumeBelow;

     mp.asurfaceAbove = p.asurfaceAbove;
     mp.asurfaceBelow = p.asurfaceBelow;

     mp.EaVolume_a = p.Ea_above;
     mp.EaVolume_b = p.Ea_below;
     //alphaOverA = alpha/Avolume
     mp.alphaVolume = p.alphaOverA * mp.AvolumeAbove;

     mp.Vso = p.Vso + asymmetry * p.VsoNZ;
     mp.Rso = std::pow(A,1./3.) * p.rso0 + p.rso1;
     mp.aso = p.aso0 + p.aso1 / std::pow(A,1./3.);
     mp.AWso = p.AWso;
     mp.BWso = p.BWso;
     mp.V_wine = p.V_wine;
     mp.R_wine = p.R_wine;
     mp.rho_wine = p.rho_wine;
         return mp;
}

NuclearParameters
read_nucleus_parameters( const std::string &filename ) {

    std::ifstream file( filename.c_str() );

    if (file.fail()) {
      std::cout << "couldn't open data file " << filename << std::endl;
      std::abort();
    }

    double Zp;
    double Z;
    double A;
    double Ef;
    int readCoulomb;

    file >> Zp >> Z >> A >> Ef >> readCoulomb;

    int Nfermi;
    file >> Nfermi;

    //FIXME this is just temporary. Should I incorporate these
    // into the NucleusParameters object?
    int Np;
    int Nh; 
    int lp;
    int lh;
    double jp;
    double jh;
    double eh;
    double ep;

    if (Nfermi < 1 || Nfermi > 2) 
        std::cout << "Nfermi not possible" << std::endl;

    file >> Nh >> jh >> lh >> eh;
    if (Nfermi == 2) file >> Np >> jp >> lp >> ep;
    else {
        Np = Nh;
        jp = jh;
        lp = lh;
        ep = eh;
    }

    double gap;
    double gapOther;

    file >> gap  >> gapOther;
    double gapMin = std::min(gap,gapOther);
    double Wgap = gap/2. + gapMin;

    return NuclearParameters( A, Z, Ef, gap, Wgap, readCoulomb, Nh, lh, jh, eh, Np, lp, jp, ep );

}

