/* Eigenvector and eigenvalue functions.
 *
 * Copyright Mark Burnett, April 2009
 */
#ifndef _NUMERIC_LINALG_EIGEN_H_
#define _NUMERIC_LINALG_EIGEN_H_

#include <utility>
#include <vector>
#include <exception>
#include "types.h"

namespace ublas = boost::numeric::ublas;

namespace numeric { namespace linalg {

// Real versions
cvector_t eigenvalues( const matrix_t &mat );
std::pair< cvector_t, matrix_t > eig( const matrix_t &mat );

// Complex versions
cvector_t eigenvalues( const cmatrix_t &mat );
std::pair< cvector_t, cmatrix_t > eig( const cmatrix_t &mat );

// Sorted eigenvalue versions
std::vector< double >
sorted_eigenvalues( const matrix_t &m );

// Complex Matrix Inversion
void 
inverse_mtx( cmatrix_t &m );

} // end namespace linalg
} // end namespace numeric

#endif // _NUMERIC_LINALG_EIGEN_H_
