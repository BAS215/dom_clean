
// Map the parameter p_external in [a, b] to a new parameters
// in [ -Inf, Inf ]
double Map_to_inf( double a, double b, double p_external ) {

    // a is the lower bound
    // b is the upper bound
    return std::asin( 2 * ( p_external - a ) / ( b - a ) - 1 );

}

// Map the parameter p_internal in [-Inf, Inf] to the original parameter 
double Map_back( double a, double b, double p_internal ) {

    // a is the lower bound
    // b is the upper bound
    return a + ( ( b - a ) / 2 ) * ( std::sin( p_internal ) + 1 );
}
