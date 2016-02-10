#ifndef _DENSITY_H_
#define _DENSITY_H_

#include <vector>
#include <cmath>
#include "numerical.h"

std::vector<double>
point_density_lj( const std::vector<double> &rmesh, const matrix_t &d_mtx_lj,
                  double xj ); 

std::vector<double>
point_density_qp_nlj( const std::vector<double> &QPF, double xj ); 

std::vector<double>
so_correction_nlj( const std::vector<double> &rmesh, 
                   const std::vector<double> &chd_nlj, 
                   int l, double xj, double tz, double occupation ); 

std::vector<double>
folded_ch_density( const std::vector<double> &rmesh, 
                   const std::vector<double> &rweights,
                   const std::vector<double> &point_dist, double tz, double A ); 
std::vector<double>
matter_distribution( const std::vector<double> &rmesh, 
                     const std::vector<double> &rweights,
                     const std::vector<double> &point_dist, 
                     double A ); 

#endif // _DENSITY_H_
