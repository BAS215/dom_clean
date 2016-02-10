// Manybody stuff
#include "manybody.h"
#include <cstdlib>
#include <fstream>


// This actually calculates r * r' * deltaR * G( r, r' )
cmatrix_t propagator( const std::vector<double> &rmesh, double E, int l, 
                      double j, const NuclearParameters &Nu, pot &U ) {

    cmatrix_t ham = c_hamiltonian( rmesh, E, l, j, Nu.A, U ); 

    // Form of Propagator is 1 / ( E - H ), where H is the hamiltonian
    ham *= -1.0;
    cmatrix_t G( ham ); // Propagator

    for( unsigned int i = 0; i < rmesh.size(); ++i ) G( i, i ) += E;
        
    numeric::linalg::inverse_mtx( G );
    return G;
}


// Find the bound states by diagonalizing Hamiltonian
// Calculates the eigenvalues only
std::vector< double >
bound_eigvals( const matrix_t &Ham ) {

    // Find eigenvalues
    cvector_t evals = numeric::linalg::eigenvalues( Ham );
    std::list< double > neg_evals;

    // Find the negative solutions
    for ( unsigned int m = 0; m < evals.size(); ++m ) {
        if ( real( evals( m ) ) < 0 ) {
            neg_evals.push_back( real( evals( m ) ) );

        }
    }
    
    // Check to make sure that eigenvalues were found
    if ( neg_evals.empty() ) {
        std::cout << "No Bound States Found!" << std::endl;
    }

    neg_evals.sort();


    // return a vector instead of a list
    std::list< double >::iterator it;
    std::vector< double > evals_vec;
    for( it = neg_evals.begin(); it != neg_evals.end(); ++it ) {

        evals_vec.push_back( *it );
    }

    return evals_vec;

}

// Find the real eigenvalues by diagonalizing Hamiltonian
// Calculates the eigenvalues only
std::vector< double >
real_eigvals( const matrix_t &Ham ) {

    // Find eigenvalues
    cvector_t evals = numeric::linalg::eigenvalues( Ham );
    std::list< double > real_evals;

    // Find the negative solutions
    for ( unsigned int m = 0; m < evals.size(); ++m ) {

        real_evals.push_back( real( evals( m ) ) );
    }
    
    // Check to make sure that eigenvalues were found
    if ( real_evals.empty() ) {
        std::cout << "No Real Eigenvalues Found!" << std::endl;
    }

    real_evals.sort();

    // return a vector instead of a list
    std::list< double >::iterator it;
    std::vector< double > evals_vec;
    for( it = real_evals.begin(); it != real_evals.end(); ++it ) {

        evals_vec.push_back( *it );
    }

    return evals_vec;

}

bool 
compare_N( std::pair< double, std::vector<double> > pair1, 
           std::pair< double, std::vector<double> > pair2 ) {

    if( pair1.first < pair2.first ) return true;
    else return false;

}

// Find the bound states by diagonalizing Hamiltonian
// Calculates the eigenvalues and the eigenvectors
std::vector< eigen_t >
bound_eigvecs( const matrix_t &Ham ) {

    // Find eigenvalues and eigenvectors
    std::pair< cvector_t, matrix_t > eval_evec = numeric::linalg::eig( Ham );
    cvector_t evals = eval_evec.first;
    matrix_t evecs = eval_evec.second;

    // find the eigenvalue corresponding to N
    std::list< eigen_t > bound_list;

    for ( unsigned int m = 0; m < evals.size(); ++m ) {
        if ( real( evals( m ) ) < 0 ) {

            std::vector< double > eigvec;
            for ( unsigned int i = 0; i < evecs.size1(); ++i ) {
                eigvec.push_back( evecs( i, m ) );
            }
            
            bound_list.push_back( 
                std::make_pair( real( evals( m ) ), eigvec ) );

        }
    }

    bound_list.sort( compare_N );

    // return a vector instead of a list
    std::list< eigen_t >::iterator it;
    std::vector< eigen_t > bound_vec;
    for ( it = bound_list.begin(); it != bound_list.end(); ++it ) {

        bound_vec.push_back( *it );
        
    }

    return bound_vec;

}

// Find the real states by diagonalizing Hamiltonian
// Calculates the eigenvalues and the eigenvectors
std::vector< eigen_t >
real_eigvecs( const matrix_t &Ham ) {

    // Find eigenvalues and eigenvectors
    std::pair< cvector_t, matrix_t > eval_evec = numeric::linalg::eig( Ham );
    cvector_t evals = eval_evec.first;
    matrix_t evecs = eval_evec.second;

    // find the eigenvalue corresponding to N
    std::list< eigen_t > real_list;

    for ( unsigned int m = 0; m < evals.size(); ++m ) {

        std::vector< double > eigvec;
        for ( unsigned int i = 0; i < evecs.size1(); ++i ) {

            eigvec.push_back( evecs( i, m ) );
        }
            
        real_list.push_back( std::make_pair( real( evals( m ) ), eigvec ) );
    }

    real_list.sort( compare_N );

    // return a vector instead of a list
    std::list< eigen_t >::iterator it;
    std::vector< eigen_t > real_vec;
    for ( it = real_list.begin(); it != real_list.end(); ++it ) {

        real_vec.push_back( *it );
        
    }

    return real_vec;

}

// Calculates the eigenvalues of the Schroedinger-like equation
// self-consistently
std::vector<eigen_t>
slfcn_eig( const std::vector<double> &rmesh, double Estart, int nb, int l, 
           double j, const NuclearParameters &Nu, pot &U ) {

    // construct Hamiltonian with starting energy value
    matrix_t ham = re_hamiltonian( rmesh, Estart, l, j, Nu.A, U );

    // get the negative eigenvalues sorted from lowest to highest
    std::vector< double > neg_evals = bound_eigvals( ham );

    // If the number of bound states for a particular lj is less
    // than the expected number, the problem may be due to a
    // level that is near the continuum. These levels can be
    // quite sensitive to the value of Estart.
    double Eth = -0.05; // Threshold energy
    while ( neg_evals.size() < nb ) {

        Estart += std::abs( Estart ) / 2;

        ham.clear();
        neg_evals.clear();

        ham = re_hamiltonian( rmesh, Estart, l, j, Nu.A, U );
        neg_evals = bound_eigvals( ham );

        if ( Estart >= Eth ) {
            std::cout << "No bound state found below E = " << Eth << std::endl;
            break;
        }
    }

    std::vector<double> Ein_vec;
    for ( unsigned int n = 0; n < neg_evals.size(); ++n ) {
        
        Ein_vec.push_back( neg_evals.at( n ) ); 
    }


    std::vector<double> Eout_vec;
    for ( unsigned int n = 0; n < Ein_vec.size(); ++n ) {
        
        ham.clear();
        neg_evals.clear();
        ham = re_hamiltonian( rmesh, Ein_vec.at( n ), l, j, Nu.A, U );
        neg_evals = bound_eigvals( ham );

        Eout_vec.push_back( neg_evals.at( n ) ); 

    }

    if ( Ein_vec.size() != Eout_vec.size() ) {

        std::cout << "# of energies in is not equal to # of energies out." 
                  << std::endl; 

        std::cout << "Ein.size() = " << Ein_vec.size() << std::endl;
        std::cout << "Eout.size() = " << Eout_vec.size() << std::endl;
        
        throw std::exception();
    }

    // find self consistent eigenvalues for each principal quantum number N
    double tol = 0.005; // tolerance
    double tol2 = 0.01;
    double dtol = 0.005; // 

    std::vector< eigen_t > eig_info; // vector for storing eigenvalues
    for ( unsigned int N = 0; N < Ein_vec.size(); ++N ) {

        double Ein = Ein_vec.at( N );
        double Eout = Eout_vec.at( N );

        int iter = 0;
        int iter_max = 30;
        int iter_med = 8;

        std::cout << "------------------------------------------" << std::endl;
        std::cout << " " << std::endl;
        std::cout << "n = " << N << " l = " << l << " j = " << j << std::endl;
        std::cout << " " << std::endl;

        while ( std::abs( Eout - Ein ) > tol ) {

            iter += 1;
            std::cout << iter << " " << Ein << " " << Eout << std::endl;

            // Assign new input energy
            if ( ( iter > iter_max ) && ( std::abs( Eout - Ein ) > tol2 ) ) {
                
                char ans;
                double new_tol;
                std::cout << "Not converging. Increase tolerance?" << std::endl;
                std::cin >> ans;
                if ( ans == 'y' ) {
                    
                    std::cout << "Please enter new tolerance..." << std::endl;
                    std::cin >> new_tol;
                    tol = new_tol;
                }
                else {

                    std::cout << "Enter Input Energy Manually:" << std::endl;
                    std::cin >> Ein;

                    if ( Ein > 0 ) std::abort();
                }
            }
            else if ( ( iter > iter_med ) && 
                      ( std::abs( Eout - Ein ) > tol2 ) ) {

                Ein = ( Eout + Ein ) / 2.;
            }
            else Ein = Eout;

            // construct Hamiltonian with new energy value
            ham.clear();
            ham = re_hamiltonian( rmesh, Ein, l, j, Nu.A, U );

            // get new negative eigenvalues
            std::vector< double > neg_evals2 = bound_eigvals( ham );
            if( neg_evals2.empty() ) continue;

            Eout = neg_evals2.at( N ); // New output energy

        } // end self-consistency loop

        ham.clear();
        ham = re_hamiltonian( rmesh, Eout, l, j, Nu.A, U );

        std::vector< eigen_t > bound_info = bound_eigvecs( ham );
        if( bound_info.empty() ) continue;

        double eval = bound_info[N].first;
        std::vector<double> &evec = bound_info[N].second;

        if ( std::abs( eval - Eout ) > tol ) {
            std::cout << "Functions 'numeric::linalg::eigenvalues' " 
                      << "and " <<  "'numeric::linalg::eig' " 
                      << "not " << "consistent with each other." << std::endl;
            
            std::cout << "From eigenvalues: " << Eout << std::endl;
            std::cout << "From eig: " << eval << std::endl;

            std::cout << "Increasing Tolerance to " << tol + dtol << std::endl;
        }

        if ( std::abs( eval - Eout ) > ( tol + dtol ) ) {
            std::cout << "Functions 'numeric::linalg::eigenvalues' " 
                      << "and " <<  "'numeric::linalg::eig' " 
                      << "not " << "consistent with each other." << std::endl;
            
            std::cout << "From eigenvalues: " << Eout << std::endl;
            std::cout << "From eig: " << eval << std::endl;

            throw std::exception();
        }

        // Normalize Eigenfunction
        std::vector<double> efx = normalize( rmesh, evec );

        eig_info.push_back( std::make_pair( eval, efx ) );
    } // end loop over N

    return eig_info;
}

// Calculates the eigenvalues of the Schroedinger-like equation
// self-consistently
double
find_level( const std::vector<double> &rmesh, double Estart, int N, int l, 
            double j, double tol, const NuclearParameters &Nu, pot &U ) {

    // construct Hamiltonian with starting energy value
    matrix_t ham = re_hamiltonian( rmesh, Estart, l, j, Nu.A, U );

    // get the negative eigenvalues sorted from lowest to highest
    std::vector< double > evals = real_eigvals( ham );

    std::vector<double> Ein_vec;
    double Ein = evals.at( N );

    ham.clear();
    evals.clear();
    ham = re_hamiltonian( rmesh, Ein, l, j, Nu.A, U );
    evals = real_eigvals( ham );

    double Eout = evals.at( N );

    // find self consistent eigenvalues for each principal quantum number N
    double tol2 = 2 * tol; 

    int iter = 0;
    int iter_max = 30;
    int iter_med = 8;

    std::cout << "------------------------------------------" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "n = " << N << " l = " << l << " j = " << j << std::endl;
    std::cout << " " << std::endl;

    while ( std::abs( Eout - Ein ) > tol ) {

        iter += 1;
        std::cout << iter << " " << Ein << " " << Eout << std::endl;

        // Assign new input energy
        if ( ( iter > iter_max ) && ( std::abs( Eout - Ein ) > tol2 ) ) {
            
            char ans;
            double new_tol;
            std::cout << "Not converging. Increase tolerance?" << std::endl;
            std::cin >> ans;
            if ( ans == 'y' ) {
                
                std::cout << "Please enter new tolerance..." << std::endl;
                std::cin >> new_tol;
                tol = new_tol;
            }
            else {

                std::cout << "Enter Input Energy Manually:" << std::endl;
                std::cin >> Ein;

                if ( Ein > 0 ) std::abort();
            }
        }
        else if ( ( iter > iter_med ) && 
                  ( std::abs( Eout - Ein ) > tol2 ) ) {

            Ein = ( Eout + Ein ) / 2.;
        }
        else Ein = Eout;

        // construct Hamiltonian with new energy value
        ham.clear();
        ham = re_hamiltonian( rmesh, Ein, l, j, Nu.A, U );

        // get new eigenvalues
        std::vector< double > evals2 = real_eigvals( ham );
        if( evals2.empty() ) continue;

        Eout = evals2.at( N ); // New output energy

    } // end self-consistency loop

    ham.clear();
    ham = re_hamiltonian( rmesh, Eout, l, j, Nu.A, U );

    std::vector< double > evals3 = real_eigvals( ham );
    double eval = evals3.at( N );

    // check stability
    if ( std::abs( eval - Eout ) > tol ) {
        std::cout << "Eigenvalue not converged within " << tol 
                  << " MeV." << std::endl;
        std::cout << "Difference is " << eval - Eout << " MeV." << std::endl;
    }

    return Eout;
}

eigen_t
find_boundstate( const std::vector<double> &rmesh, double Estart, int N, int l, 
                 double j, double tol, const NuclearParameters &Nu, pot &U ) {

        // get level
        double QPE = find_level( rmesh, Estart, N, l, j, tol, Nu, U );
        matrix_t ham = re_hamiltonian( rmesh, QPE, l, j, Nu.A, U );

        std::vector< eigen_t > eig_info = real_eigvecs( ham );

        double eval = eig_info[N].first;
        std::vector<double> &evec = eig_info[N].second;

        if ( std::abs( eval - QPE ) > tol ) {
            std::cout << "Functions 'numeric::linalg::eigenvalues' " 
                      << "and " <<  "'numeric::linalg::eig' " 
                      << "not " << "consistent with each other." << std::endl;
            
            std::cout << "From eigenvalues: " << QPE << std::endl;
            std::cout << "From eig: " << eval << std::endl;

            // choose level that is outside continuum
            ham.clear();
            ham = re_hamiltonian( rmesh, eval, l, j, Nu.A, U );
            std::vector< eigen_t > eig_info2 = real_eigvecs( ham );
            double eval2 = eig_info2[N].first;
            std::vector< double > &evec2 = eig_info2[N].second;

            if ( eval2 > eval ) {

                for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

                    evec[i] = evec2[i];
                }
                eval = eval2;
            }

            std::cout << "New Calculation: " << eval << std::endl;

        }

        // Normalize Eigenfunction
        std::vector<double> efx = normalize( rmesh, evec );

        return std::make_pair( eval, efx );
}

std::vector< double > 
normalize( const std::vector< double > &rmesh, 
           const std::vector<double> &eigenvector ) {

    double delt = rmesh.at(1) - rmesh.at(0);

    double norm2 = 0;
    for( unsigned int i = 0; i < rmesh.size(); ++i ) {

        norm2 += std::pow( eigenvector[ i ], 2 ) * delt;
        
    }
    double norm = std::sqrt( norm2 );

    std::vector< double > efx;
    double sign = 1;
    if ( eigenvector[0] < 0 ) sign = -1;
    for( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        efx.push_back( sign * eigenvector[i] / norm / rmesh[i] );
    }

    return efx;
}


// Calculate the Spectroscopic Factor. 'QPE' is the quasiparticle
// ( or quasihole ) energy and 'QPF' is the corresponding overlap function
double sfactor( const std::vector<double> &rmesh, double QPE, int l, 
                double xj, const std::vector<double> &QPF, pot &U ) {

    U.setEnergy( QPE );
    U.setAM( l, xj );

    double delt = rmesh.at(1) - rmesh.at(0);

    double dsigma_local = 0;
    double dsigma_nonlocal = 0;

    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        
        // local parts
        double dV_local = U.der_disp_localPart( rmesh[i] );

        dsigma_local += delt * std::pow( rmesh[i], 2 ) * dV_local 
                      * std::pow( QPF[i], 2 );
                      
        // nonlocal parts
        for ( unsigned int k = 0; k < rmesh.size(); ++k ) {
            
            double dV_nonlocal = 
                U.der_disp_nonlocalPart( rmesh[i], rmesh[k] );

            // nonlocal part is already multiplied by rr'
            dsigma_nonlocal += delt * delt * dV_nonlocal
                             * rmesh[i] * QPF[i] * rmesh[k] * QPF[k];
        }
    }

    std::cout << "dsigma_local = " << dsigma_local << std::endl;
    std::cout << "dsigma_nonlocal = " << dsigma_nonlocal << std::endl;

    double dsigma = dsigma_local + dsigma_nonlocal;
    return 1.0 / ( 1.0 - dsigma );
}

double rms_radius( const std::vector<double> &rmesh, 
                   const std::vector<double> &QPF ) {

    double delt = rmesh.at(1) - rmesh.at(0);

    double rsq = 0;

    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {

        rsq += std::pow( rmesh[i], 4 ) * delt * std::pow( QPF[i], 2 );
    }

    return std::sqrt( rsq );
}

// Occupation calculated from the density matrix. 'd_mtx' is the 
// density matrix
double occupation( const std::vector<double> &rmesh, const matrix_t &d_mtx, 
                   const std::vector<double> &QPF ) {


    double delt = rmesh.at(1) - rmesh.at(0);

    double occ = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        for( unsigned int j = 0; j < rmesh.size(); ++j ) {

            occ += std::pow( rmesh[i], 2 ) * QPF[i] * d_mtx( i, j )
                 * std::pow( rmesh[j], 2 ) * QPF[j] * std::pow( delt, 2 );
        }
    }

    return occ;
}

double 
width( const std::vector<double> &rmesh, double QPE, int L, double J,
       const std::vector< double > &QPF, pot &U ) {

    U.setEnergy( QPE );
    U.setAM( L, J );

    double delt = rmesh.at( 1 ) - rmesh.at( 0 );

    double sum = 0;
    for ( unsigned int i = 0; i < rmesh.size(); ++i ) {
        for ( unsigned int j = 0; j < rmesh.size(); ++j ) {

            // potential is already multiplied by rr'
            sum += delt * delt * rmesh[i] * rmesh[j] 
                 * QPF[i] * QPF[j] * U.nonlocalIM( rmesh[i], rmesh[j] );
        }
    }

    return 2 * std::abs( sum );
}

