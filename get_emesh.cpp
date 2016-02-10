
std::vector< mesh_t >
boundRspace::get_emeshes( const std::vector< double > &rmesh, 
                          double Emin, double Emax, 
                          const std::vector< lj_eigen_t > &bound_levels ) {

    std::cout << "Creating energy meshes ... " << std::endl;

    // vector for storing the energy meshes for each LJ combination
    std::vector< mesh_t > emesh_vec;

    // loop over L
    for ( int L = 0; L < Lmax + 1; ++L ) {

        // loop over J
        for ( int up = -1; up < 2; up +=2 ) {

            double J = L + up / 2.0;
            if ( J < 0 ) continue;

            int lj_index = index_from_LJ( L, J );
            const std::vector< eigen_t > &levels = bound_levels.at( lj_index );

            lj_eigen_t qh_in_continuum_vec;
            for ( unsigned int N = 0; N < levels.size(); ++N ) {

                if ( levels[N].first >= 0 ) { 
                    
                    std::cout << "Error: level >= 0." << std::endl;
                    std::cout << "N, L, J = " << N << " " << L << " " << J
                              << std::endl;

                    std::cout << "Level = " << levels[N].first << std::endl;
                    
                    std::abort();
                }

                // the locations of the quasihole peaks that
                // are in or near the continuum will be used
                // to construct the energy mesh. 
                if ( levels[N].first < ( Emax + 3 ) ) {
                    
                    qh_in_continuum_vec.push_back( levels[N] );

                    if ( levels[N].first < Emax ) {
//
//                        std::cout << " " << std::endl;
//                        std::cout << "Level in continuum." << std::endl;
//                        std::cout << " " << std::endl;
                    }
                    else {

//                        std::cout << " " << std::endl;
//                        std::cout << "Level outside but near continuum." 
//                                 << std::endl;
//                       std::cout << " " << std::endl;
                    }
                }

            } // End loop over N

            // Create Energy Grid for calculating density matrix etc.
            // based on the location of the quasihole peaks
            double Elowmax = -60; // default
            double Emedmax = Emax - ph_gap; // default
            double Edeltlow = 5.; // default
            double Edeltmed = 1;// default
            double Edelthigh = 0.2; // default

            if ( qh_in_continuum_vec.size() == 1 ) {

                double QPE = qh_in_continuum_vec[0].first;
                std::vector< double > &QPF = qh_in_continuum_vec[0].second;

                // get rough estimate of quasihole width
                double gamma = approx_width( rmesh, QPE, L, J, QPF );

                if ( QPE > Emedmax ) {
                    
                    Edelthigh = 0.05;
                }
                else {

                    Elowmax = QPE - gamma;
                    Edeltmed = 0.2;
                }

            }
            else if ( qh_in_continuum_vec.size() > 1 ) {

                // Find the minimum and maximum bound-state
                // energies located in or near the continuum
                eigen_t &bound_max = qh_in_continuum_vec.back();
                eigen_t &bound_min = qh_in_continuum_vec.front();

                double QPE_max = bound_max.first;
                double QPE_min = bound_min.first;

                if ( QPE_max > Emedmax ) Edelthigh = 0.05;

                // get rough estimate of quasihole width of 
                // lowest level
                double gamma = 
                    approx_width( rmesh, QPE_min, L, J, bound_min.second );

                Elowmax = QPE_min - gamma;
                Edeltmed = 0.2 ;

            }

            mesh_t emesh = energy_mesh_fast( Emin, Elowmax, Emedmax, Emax, 
                                             Edeltlow, Edeltmed, Edelthigh );

            emesh_vec.push_back( emesh );

        } // End loop over J
    } // End loop over L

    if ( emesh_vec.size() != static_cast< unsigned >( 2 * Lmax + 1 ) ) {

        std::cout << "Some  LJ combinations are missing." << std::endl;
        std::cout << "emesh_vec.size() = " << emesh_vec.size() << std::endl;
        std::cout << "expected size = " << 2 * Lmax + 1 << std::endl;
    }

    std::cout << "Finished creating energy meshes." << std::endl;

    return emesh_vec;
}
