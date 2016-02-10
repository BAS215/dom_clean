//
#include "meshes.h"

/*
std::vector<double>
energy_mesh( const std::vector<double> &points, 
             const std::vector<double> &widths ) {

    if ( points.size() < 2 ) {
        
        std::cout << "Must have at least two points to define mesh" 
                  << std::endl;

        throw std::exception();
    }

    for ( unsigned int p = 0; p < points.size() - 1; ++p ) {

        double low = points[p];
        double high = points[p+1];

        int epts = static_cast<int>( ( high - low ) / widths[p] );
    }
}
*/

std::vector<std::pair< double, double > >
energy_mesh_slow( double Emin, double Emax, double Edelt ) {

    int epts = static_cast<int>( ( Emax - Emin ) / Edelt );

    std::vector< std::pair< double, double > > emesh;
    for ( int i = 0; i < epts; ++i ) {

        double E = Emin + i * Edelt;

        emesh.push_back( std::make_pair( E, Edelt ) );
    }

    return emesh;
}

std::vector< std::pair< double, double > >
energy_mesh_fast( double Emin, double Elowmax, double Emedmax, double Emax,
                  double Edeltlow, double Edeltmed, double Edelthigh ) {

    std::vector<std::pair< double, double > > emesh;

    // Very Negative Energy
    int epts = static_cast<int>( ( Elowmax - Emin ) / Edeltlow );

    for ( int i = 0; i < epts; ++i ) {

        double E = Emin + ( i + 0.5 ) * Edeltlow;
        emesh.push_back( std::make_pair( E, Edeltlow ) ); 
    }

    // Medium Negative Energy
    double Emedlow = emesh.back().first + Edeltlow / 2;
    epts = static_cast<int>( ( Emedmax - Emedlow ) / Edeltmed );

    for ( int i = 0; i < epts; ++i ) {

        double E = Emedlow + ( i + 0.5 ) * Edeltmed;
        emesh.push_back( std::make_pair( E, Edeltmed ) );
    }

    // Energy Near Cutoff
    double Ecutlow = emesh.back().first + Edeltmed / 2;

    //change Ecutlow until ( Emax - Ecutlow ) / Edelthigh is an integer
    epts = static_cast<int>( ( Emax - Ecutlow ) / Edelthigh );

    double Ecutlow_new = Emax - epts * Edelthigh;

    for ( int i = 0; i < epts; ++i ) {
        
        double E = Ecutlow_new + ( i + 0.5 ) * Edelthigh;
        emesh.push_back( std::make_pair( E, Edelthigh ) );
    }

    //std::cout << "E_last = " << emesh.back().first << std::endl;
    //std::cout << "E_max = " << Emax << std::endl;

    return emesh;
}
