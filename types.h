// Copyright 2009 Mark Burnett
// Numeric types used in the code will be defined here.

#ifndef _NUMERIC_TYPES_H_
#define _NUMERIC_TYPES_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <vector>
#include <complex>

namespace ublas = boost::numeric::ublas;

// Real value types
typedef ublas::vector< double > vector_t;
typedef ublas::matrix< double, ublas::column_major > matrix_t;

// Complex valued types
typedef std::complex< double > complex_t;
typedef ublas::vector< complex_t > cvector_t;
typedef ublas::matrix< complex_t, ublas::column_major > cmatrix_t;

#endif // _NUMERIC_TYPES_H_
