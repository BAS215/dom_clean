//
#ifndef _MESHES_H_
#define _MESHES_H_

#include <utility>
#include <vector>
#include <iostream>

std::vector<std::pair< double, double > >
energy_mesh_slow( double Emin, double Emax, double Edelt );

std::vector< std::pair< double, double > >
energy_mesh_fast( double Emin, double Elowmax, double Emedmax, double Emax,
                  double Edeltlow, double Edeltmed, double Edelthigh );

#endif //_MESHES_H_
