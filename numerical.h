//

#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_


#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <exception>
#include "types.h"

double der_2pt( double delta, double f0, double f_plus ); 
double der_3pt( double delta, double f_minus, double f_plus ); 
double der_4pt( double delta, double f_minus2, double f_minus,
                double f_plus, double f_plus2 ); 

std::pair< double, double > 
lin_regression( const std::vector<double> &x_vec,
                const std::vector<double> &y_vec );

#endif // _NUMERICAL_H_
