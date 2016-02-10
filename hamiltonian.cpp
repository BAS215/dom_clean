//
#include "hamiltonian.h"

// Wrapper for Bob's potential class.
// This function returns a potential class with 
// all the appropriate initializations.
// >> type indicates form factor used for nonlocal potentials:
//      * 1 is the form of Perey and Buck (average)
//      * 0 is the form of D. Van Neck
// >> mvolume is used in the form for the energy dependence 
//    of the imaginary volume potential. 
// >> AsyVolume is used to specify whether the imaginary volume
// >> potential will have an asymmetry dependence (1 == yes, 0 == no)
// >> tz is the isospin of projectile (-0.5 for neutrons, 0.5 for protons)
// >> Nu holds information about the target (Fermi energy, A, Z, etc.)
// >> p is a struct holding all the parameters
pot get_bobs_pot( int type, int mvolume, int AsyVolume, double tz, 
		  const NuclearParameters &Nu, const Parameters &p ) {

    double Zp; 
    if ( tz > 0 ) Zp = 1;
    else Zp = 0;

    pot Pot( type );
    Pot.init( Nu.Z, Zp, Nu.A, Nu.readCoulomb );
    Pot.load( p.Rc, p.VHFvol, p.VHFsur, p.RHF, p.aHF,  p.RHFs, p.aHFs, p.beta_nl_R0, p.AHF,
              p.beta_nl_R1, p.RsurfaceAbove,p.RsurfaceBelow, p.asurfaceAbove,p.asurfaceBelow,
	      p.AsurfaceAbove, p.AsurfaceBelow, p.BsurfaceA, p.CsurfaceA, p.DsurfaceA, p.Bsurface, p.Csurface, p.Dsurface, 
              Nu.Wgap * p.fGap_A , Nu.Wgap * p.fGap_B , Nu.Ef, p.beta_nl_I0,
              p.beta_nl_I1,  p.beta_nl_I0_sur,  p.beta_nl_I1_sur, 
	      p.RvolumeAbove,p.RvolumeBelow, p.deltaRvolume, p.expRvolume, p.avolumeAbove,p.avolumeBelow,
              p.AvolumeAbove, p.AvolumeBelow, p.BvolumeAbove,p.BvolumeBelow,
	      p.EpvolumeAbove,p.EpvolumeBelow, mvolume,
              AsyVolume, p.alphaVolume, p.EaVolume_a , p.EaVolume_b, p.Rso,
              p.aso, p.Vso, p.AWso, p.BWso, p.V_wine, p.R_wine, p.rho_wine );


    return Pot;
}

// returns the real part of the Hamiltonian in coordinate space.
// This is needed when solving the Schroedinger-like equation 
// for the overlap functions.
//
matrix_t 
re_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J, double A, pot &U ) {

    U.setEnergy( Ecm ); 
    U.setAM( L, J );

    double rdelt = rmesh.at(1) - rmesh.at( 0 );

    // initialize matrix
    matrix_t ham( rmesh.size(), rmesh.size() );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    for ( unsigned int j = 0; j <= i; ++j ) {
        
        ham( i, j ) = 0.0;
        ham( j, i ) = ham( i, j );
    }
    }

    // \hbar^2 / ( 2 * mu ), mu is reduced mass
    double mu_factor = A / ( A + 1.0 );
    double fac = - 1 / ( U.kconstant * mu_factor );

    // Construct Tridiagonal matrix for kinetic energy term
    // and add in diagonal parts of potential
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        // Diagonal Part of Kinetic Energy
        double KE_diagonal = -2 * fac / std::pow( rdelt, 2 )
                          - fac * L * ( L + 1 ) / std::pow( rmesh[i], 2 );

        // Off-diagonal part of Kinetic energy
        double KE_off_diagonal = fac / std::pow( rdelt, 2 );

        // Diagonal Parts of Potential Energy
        double PE_diagonal = real ( U.localPart( rmesh[i] ) );

        ham( i, i ) = PE_diagonal + KE_diagonal;
        if ( i != 0 ) ham( i, i - 1 ) = KE_off_diagonal;
        if ( i != rmesh.size() - 1 ) ham( i, i + 1 ) = KE_off_diagonal;

        // Boundary Condition
        if ( i == 0 ) ham( 0, 0 ) -= KE_off_diagonal * std::pow( -1.0, L );

        // Now add in nonlocal components
        for ( unsigned int j = 0; j <= i; ++j ) {

            // nonlocalPart is already multiplied by r, r'
            ham( i, j ) += real( U.nonlocalPart( rmesh[i], rmesh[j] ) ) * rdelt;

            if ( j != i ) ham( j, i ) = ham( i, j );

            //std::cout << i << " " << j << " " << ham( i, j ) << std::endl;
        } // end loop over j

    } // end loop over i

    return ham;
}

// returns the complex part of the Hamiltonian. This is needed when
// solving the Dyson equation to get the propagator 
cmatrix_t 
c_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J, double A, pot &U ) {

    U.setEnergy( Ecm ); 
    U.setAM( L, J );

    double rdelt = rmesh.at(1) - rmesh.at( 0 );

    // initialize matrix
    cmatrix_t ham( rmesh.size(), rmesh.size() );
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
    for ( unsigned int j = 0; j <= i; ++j ) {
        
        ham( i, j ) = complex<double>( 0.0, 0.0 );
        ham( j, i ) = ham( i, j );
    }
    }

    // \hbar^2 / ( 2 * mu ), mu is reduced mass
    double mu_factor = A / ( A + 1.0 );
    double fac = - 1 / ( U.kconstant * mu_factor );

    // Construct Tridiagonal matrix for kinetic energy term
    // and add in diagonal parts of potential
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        // Diagonal Part of Kinetic Energy
        double KE_diagonal = -2 * fac / std::pow( rdelt, 2 )
                          - fac * L * ( L + 1 ) / std::pow( rmesh[i], 2 );

        // Off-diagonal part of Kinetic energy
        double KE_off_diagonal = fac / std::pow( rdelt, 2 );

        // Diagonal Parts of Potential Energy
        complex<double> PE_diagonal = U.localPart( rmesh[i] );

        ham( i, i ) = PE_diagonal + KE_diagonal;
        if ( i != 0 ) ham( i, i - 1 ) = KE_off_diagonal;
        if ( i != rmesh.size() - 1 ) ham( i, i + 1 ) = KE_off_diagonal;

        // Boundary Condition
        if ( i == 0 ) ham( 0, 0 ) -= KE_off_diagonal * std::pow( -1.0, L );

        // Now add in nonlocal components
        for ( unsigned int j = 0; j <= i; ++j ) {

            // nonlocalPart is already multiplied by r, r'
            ham( i, j ) += U.nonlocalPart( rmesh[i], rmesh[j] ) * rdelt;

            if ( j != i ) ham( j, i ) = ham( i, j );

        } // end loop over j

    } // end loop over i

    return ham;
}
