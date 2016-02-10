//
#include "numerical.h"

/* Derivative Formulas */


// 2 point formula
double der_2pt( double delta, double f0, double f_plus ) {

    return ( f_plus - f0 ) / delta;
}

// 3 point formula
double der_3pt( double delta, double f_minus, double f_plus ) {

    return ( f_plus - f_minus ) / ( 2 * delta );
}

// 4 point formula 
double der_4pt( double delta, double f_minus2, double f_minus,
                double f_plus, double f_plus2 ) {

    return ( f_minus2 - 8 * f_minus + 8 * f_plus - f_plus2 ) / ( 12 * delta );
}

// Simple Regression Model
std::pair< double, double >
lin_regression( const std::vector<double> &x_vec, 
                const std::vector<double> &y_vec ) {

    if( x_vec.size() != y_vec.size() ) {
        std::cout << "y vector and x vector need to have the same size " 
                  << "in lin_regression" << std::endl;

        throw std::exception();
    }

    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumx2 = 0;
    unsigned int N = x_vec.size();

    for( unsigned int i = 0; i < N; ++i ) {
        
        sumx += x_vec[i];
        sumy += y_vec[i];
        sumxy += x_vec[i] * y_vec[i];
        sumx2 += x_vec[i] * x_vec[i];

    }
    double denominator = N * sumx2 - std::pow( sumx, 2 );
    double m = ( sumy * sumx2 - sumx * sumxy ) / denominator; // intercept
    double b = ( N * sumxy - sumx * sumy ) / denominator; // slope

    return std::make_pair( m, b );
}
