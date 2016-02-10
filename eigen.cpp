/* Eigenvalue/vector wrappers.
 * 
 * Copyright Mark Burnett, November 2008
 */

#include "eigen.h"

#include <boost/foreach.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

namespace lapack = boost::numeric::bindings::lapack;
namespace atlas = boost::numeric::bindings::atlas;

namespace numeric { namespace linalg {

// --------------------------------------------------------------------
// Linear algebra stuff
// --------------------------------------------------------------------
// Real valued solvers
cvector_t eigenvalues( const matrix_t &mat ) {
    cvector_t vals(  mat.size1() );
//    matrix_t right_vecs( mat.size1(), mat.size2() );
    matrix_t * dummy = 0;
    matrix_t temp(mat);
    lapack::geev( temp, vals, dummy, dummy, lapack::optimal_workspace() );
    return vals; }

std::pair< cvector_t, matrix_t > eig( const matrix_t &mat ) {
    cvector_t vals(  mat.size1() );
    matrix_t right_vecs( mat.size1(), mat.size2() );
    matrix_t * dummy = 0;
    matrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return std::make_pair( vals, right_vecs ); }

// Complex valued solvers
cvector_t eigenvalues( const cmatrix_t &mat ) {
    cvector_t vals(  mat.size1() );
//    cmatrix_t right_vecs( mat.size1(), mat.size2() );
    cmatrix_t * dummy = 0;
    cmatrix_t temp(mat);
    lapack::geev( temp, vals, dummy, dummy, lapack::optimal_workspace() );
    return vals; }

std::pair< cvector_t, cmatrix_t > eig( const cmatrix_t &mat ) {
    cvector_t vals(  mat.size1() );
    cmatrix_t right_vecs( mat.size1(), mat.size2() );
    cmatrix_t * dummy = 0;
    cmatrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return std::make_pair( vals, right_vecs ); }

// Returns the real parts of the eigenvalues, sorted.
std::vector< double >
sorted_eigenvalues( const matrix_t &m ) {
    cvector_t vals = eigenvalues( m );
    std::vector< double > results;
    BOOST_FOREACH( const complex_t &v, vals ) {
        results.push_back( v.real() ); }
    std::sort( results.begin(), results.end() );
    return results; }

/*
// Returns both eigenvalues and eigenvectors, sorted by the eigenvalues.
std::vector< std::pair< complex_t, cvector_t > >
sorted_eigenvectors( const matrix_t &m ) {
    // Perform eigenvalue solution
    std::pair< cvector_t, cmatrix_t > full_pair = eig( m );
    cvector_t &vals = full_pair.first;
    cmatrix_t &vecs = full_pair.second;

    // Reorganize the solutions
    std::vector< std::pair< complex_t, cvector_t > >
        results;
    for ( int i = 0; i < boost::numeric_cast<int>(vals.size()); ++i ) {
        cvector_t v = ublas::matrix_column< cmatrix_t >( vecs, i );
        results.push_back( std::make_pair( vals(i), v ) ); }

    // Sort and returns
    std::sort( results.begin(), results.end() );
    return results; }
*/

// Returns the inverse of a complex matrix
void 
inverse_mtx( cmatrix_t &m ) {

    std::vector<int> ipiv ( m.size1() );
    int ierr = atlas::getrf( m, ipiv );
    if ( ierr == 0 ) {
        atlas::getri( m, ipiv );
    }
    else {
        std::cout << "Matrix is singular " << std::endl;
        throw std::exception();
    }
}

} // end namespace linalg
} // end namespace numeric
