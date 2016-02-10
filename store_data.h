#ifndef _STORE_DATA_
#define _STORE_DATA_

#include "boost/tuple/tuple.hpp"
#include "types.h"
#include <string>
#include <vector>
#include <utility>

struct Data_set {

    Data_set( const std::vector<double> &x_values0,
              const std::vector<double> &y_values0,
              const std::vector<double> &errors0 ) 
        : x_values( x_values0 ), y_values( y_values0 ), errors( errors0 ) {}


    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> errors;
};

struct Data_set2 {

    Data_set2( const std::vector<double> &mesh0, const matrix_t &mtx0,
               const std::vector<double> &errors0 ) 
        : mesh( mesh0 ), mtx( mtx0 ), errors( errors0 ) {}


    std::vector<double> mesh;
    matrix_t mtx;
    std::vector<double> errors;
};

typedef std::vector<Data_set> data_sets_t;
typedef std::pair< double, Data_set > data_set_pair;
typedef std::vector< data_set_pair > data_sets_2_t;
typedef std::vector<Data_set2> data_sets_mtx_t;

struct Data {

    Data( double A0, double Z0, double tz0, double Ef0,
          const Data_set &Int_form_set0,
          const Data_set &Diag_form_set0 ) 
        : A( A0 ), Z( Z0 ), tz( tz0 ), Ef( Ef0 ), 
          Int_form_set( Int_form_set0 ),
          Diag_form_set( Diag_form_set0 ) {}

    // Constants
    double A;
    double Z;
    double tz;
    double Ef;

    // Data containers
    Data_set Int_form_set;
    Data_set Diag_form_set;

};

typedef std::vector< Data > Data_vec;

#endif // _STORE_DATA_

