#ifndef _NUCLEAR_PARAMETERS_H_
#define _NUCLEAR_PARAMETERS_H_

#include <cmath>

struct NuclearParameters {

    NuclearParameters() {}
    NuclearParameters( double nA, double nZ, double nEf, double nph_gap, 
                       double nWgap, int nreadCoulomb, int nValenceHoleN, 
                       int nValenceHoleL, double nValenceHoleJ, 
                       double nValenceHoleEn, int nValenceParticleN, 
                       int nValenceParticleL, double nValenceParticleJ, 
                       double nValenceParticleEn )  
                : A( nA ), Z( nZ ), Ef( nEf ), ph_gap( nph_gap ),
                  Wgap( nWgap ), readCoulomb( nreadCoulomb ), 
                  ValenceHoleN( nValenceHoleN ),
                  ValenceHoleL( nValenceHoleL ),
                  ValenceHoleJ( nValenceHoleJ ), 
                  ValenceHoleEn( nValenceHoleEn ),
                  ValenceParticleN( nValenceParticleN ),
                  ValenceParticleL( nValenceParticleL ),
                  ValenceParticleJ( nValenceParticleJ ),
                  ValenceParticleEn( nValenceParticleEn ) { }

    double A; // number of nucleons in target
    double Z; // number of protons in target
    double Ef; // Fermi Energy
    double ph_gap; // particle-hole gap
    double Wgap; // gap used to determine where imaginary part starts
    int readCoulomb; // integer to specify whether or not to read in 
                     // exp. charge density

    int ValenceHoleN;
    int ValenceHoleL;
    double ValenceHoleJ;
    double ValenceHoleEn;

    int ValenceParticleN;
    int ValenceParticleL;
    double ValenceParticleJ;
    double ValenceParticleEn;

    double N() const { return A - Z; } // number of neutrons in target
    double Tz() const { return 0.5 * ( 2*Z - A ); }
    double T() const { return std::abs( Tz() ); }

};

#endif // _NUCLEAR_PARAMETERS_H_
