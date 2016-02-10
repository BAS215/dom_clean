//
#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "read_parameters.h"
#include "types.h"
#include "pot.h"

pot get_bobs_pot( int type, int mvolume, int AsyVolume, double tz,
                  const NuclearParameters &Nu, const Parameters &p ); 
matrix_t 
re_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J, double A, pot &U ); 
cmatrix_t 
c_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J, double A, pot &U ); 

#endif // _HAMILTONIAN_H_
