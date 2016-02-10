// Get the experimental charge density
#include <boost/math/special_functions/bessel.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>

int main() {

    std::string nucleus_string;
    int A;
    int Z;
    std::string param_string;

    std::cout << "Enter element symbol (lowercase) and A:" << std::endl;
    std::cin >> nucleus_string;

    std::cout << "Enter A:" << std::endl;
    std::cin >> A;

    std::cout << "Enter Z:" << std::endl;
    std::cin >> Z;

    std::cout << "Indicate sum-of-Guassians (sog) parametrization " <<
                 "or Fourier-Bessel (fb) parametrization:" << std::endl;
    std::cin >> param_string;

    //std::string sog_string = "sog";
    //std::string fb_string = "fb";

    /* Parameters for using the Fourier-Bessel method */
    // Numbers are from Atomic Data and Nuclear Data Tables 36, 495-536 (1987)

    // For Fourier Bessel Analysis
    std::vector<double> a_coef;
    double R_cutoff;
    double delta_Rrms; // uncertainty of root mean square radius in fm

    // For Sum of Gaussians Analysis
    std::vector<double> R_coef;
    std::vector<double> Q_coef;
    std::vector<double> A_coef;
    double gamma = 1; // for Tin isotopes (see ref. mentioned above)
    double error = 0.01;

    if( ( A == 12 ) && ( Z == 6 ) && param_string.compare("fb") == 0 ) {

        a_coef.push_back(  0.15721e-1 ); // a1
        a_coef.push_back(  0.38732e-1 ); // a2
        a_coef.push_back(  0.36808e-1 ); // a3 
        a_coef.push_back(  0.14671e-1 ); // a4
        a_coef.push_back( -0.43277e-2 ); // a5

        a_coef.push_back( -0.97752e-2 ); // a6
        a_coef.push_back( -0.68908e-2 ); // a7
        a_coef.push_back( -0.27631e-2 ); // a8
        a_coef.push_back( -0.63568e-3 ); // a9
        a_coef.push_back(  0.71809e-4 ); // a10

        a_coef.push_back(  0.18441e-3 ); // a11
        a_coef.push_back(  0.75066e-4 ); // a12
        a_coef.push_back(  0.51069e-4 ); // a13
        a_coef.push_back(  0.14308e-4 ); // a14
        a_coef.push_back(  0.23170e-5 ); // a15

        a_coef.push_back(  0.68465e-6 ); // a16

        R_cutoff = 8.0;

        delta_Rrms = 0.01; 

        error = 0.01;
    }
    else if( ( A == 12 ) && ( Z == 6 ) && param_string.compare("sog") == 0 ) {

        R_coef.push_back( 0.0 ); // R1
        R_coef.push_back( 0.4 ); // R2
        R_coef.push_back( 1.0 ); // R3
        R_coef.push_back( 1.3 ); // R4
        R_coef.push_back( 1.7 ); // R5
        R_coef.push_back( 2.3 ); // R6
        R_coef.push_back( 2.7 ); // R7
        R_coef.push_back( 3.5 ); // R8
        R_coef.push_back( 4.3 ); // R9
        R_coef.push_back( 5.4 ); // R10
        R_coef.push_back( 6.7 ); // R11

        Q_coef.push_back( 0.016690 ); // Q1
        Q_coef.push_back( 0.050325 ); // Q2
        Q_coef.push_back( 0.128621 ); // Q3
        Q_coef.push_back( 0.180515 ); // Q4
        Q_coef.push_back( 0.219097 ); // Q5
        Q_coef.push_back( 0.278416 ); // Q6
        Q_coef.push_back( 0.058779 ); // Q7
        Q_coef.push_back( 0.057817 ); // Q8
        Q_coef.push_back( 0.007739 ); // Q9
        Q_coef.push_back( 0.002001 ); // Q10
        Q_coef.push_back( 0.000007 ); // Q11

        gamma = std::sqrt( 2 / 3.0 ) * 1.20;
        error = 0.01;
        delta_Rrms = 0.003;
        R_cutoff = 9.0;

        for ( unsigned int i = 0; i < R_coef.size(); ++i ) {

            double An = Z * Q_coef[i] 
                          / ( 2 * std::pow( M_PI, 1.5 ) * std::pow( gamma, 3 ) 
                            * ( 1 + 2 * std::pow( R_coef[i] / gamma, 2 ) ) );

            A_coef.push_back( An );
        }
    }

    else if( ( A == 39 ) && ( Z == 19 ) && param_string.compare("sog") == 0 ) {

        R_coef.push_back( 0.4 ); // R1
        R_coef.push_back( 0.9 ); // R2
        R_coef.push_back( 1.7 ); // R3
        R_coef.push_back( 2.1 ); // R4
        R_coef.push_back( 2.6 ); // R5
        R_coef.push_back( 3.2 ); // R6
        R_coef.push_back( 3.7 ); // R7
        R_coef.push_back( 4.2 ); // R8
        R_coef.push_back( 4.7 ); // R9
        R_coef.push_back( 5.5 ); // R10
        R_coef.push_back( 5.9 ); // R11
        R_coef.push_back( 6.9 ); // R12

        Q_coef.push_back( 0.043308 ); // Q1
        Q_coef.push_back( 0.036283 ); // Q2
        Q_coef.push_back( 0.110517 ); // Q3
        Q_coef.push_back( 0.147676 ); // Q4
        Q_coef.push_back( 0.189541 ); // Q5
        Q_coef.push_back( 0.274173 ); // Q6
        Q_coef.push_back( 0.117691 ); // Q7
        Q_coef.push_back( 0.058273 ); // Q8
        Q_coef.push_back( 0.000006 ); // Q9
        Q_coef.push_back( 0.021380 ); // Q10
        Q_coef.push_back( 0.000002 ); // Q11
        Q_coef.push_back( 0.001145 ); // Q12

        gamma = std::sqrt( 2 / 3.0 ) * 1.45;
        error = 0.01;
        delta_Rrms = 0.003;
        R_cutoff = 10.0;

        for ( unsigned int i = 0; i < R_coef.size(); ++i ) {

            double An = Z * Q_coef[i] 
                          / ( 2 * std::pow( M_PI, 1.5 ) * std::pow( gamma, 3 ) 
                            * ( 1 + 2 * std::pow( R_coef[i] / gamma, 2 ) ) );

            A_coef.push_back( An );
        }
    }        

    else if( ( A == 40 ) && ( Z == 20 ) && param_string.compare("sog") == 0 ) {

        R_coef.push_back( 0.4 ); // R1
        R_coef.push_back( 1.2 ); // R2
        R_coef.push_back( 1.8 ); // R3
        R_coef.push_back( 2.7 ); // R4
        R_coef.push_back( 3.2 ); // R5
        R_coef.push_back( 3.6 ); // R6
        R_coef.push_back( 4.3 ); // R7
        R_coef.push_back( 4.6 ); // R8
        R_coef.push_back( 5.4 ); // R9
        R_coef.push_back( 6.3 ); // R10
        R_coef.push_back( 6.6 ); // R11
        R_coef.push_back( 8.1 ); // R12

        Q_coef.push_back( 0.042870 ); // Q1
        Q_coef.push_back( 0.056020 ); // Q2
        Q_coef.push_back( 0.167853 ); // Q3
        Q_coef.push_back( 0.317962 ); // Q4
        Q_coef.push_back( 0.155450 ); // Q5
        Q_coef.push_back( 0.161897 ); // Q6
        Q_coef.push_back( 0.053763 ); // Q7
        Q_coef.push_back( 0.032612 ); // Q8
        Q_coef.push_back( 0.004803 ); // Q9
        Q_coef.push_back( 0.004541 ); // Q10
        Q_coef.push_back( 0.000015 ); // Q11
        Q_coef.push_back( 0.002218 ); // Q12

        gamma = std::sqrt( 2 / 3.0 ) * 1.45;
        error = 0.01;
        delta_Rrms = 0.003;
        R_cutoff = 10.0;

        for ( unsigned int i = 0; i < R_coef.size(); ++i ) {

            double An = Z * Q_coef[i] 
                          / ( 2 * std::pow( M_PI, 1.5 ) * std::pow( gamma, 3 ) 
                            * ( 1 + 2 * std::pow( R_coef[i] / gamma, 2 ) ) );

            A_coef.push_back( An );
        }
    }        

    // Calcium 40
    else if( ( A == 40 ) && ( Z == 20 ) ) {
        a_coef.push_back(  0.44846e-1 ); // a1
        a_coef.push_back(  0.61326e-1 ); // a2
        a_coef.push_back( -0.16818e-2 ); // a3 
        a_coef.push_back( -0.26217e-1 ); // a4
        a_coef.push_back( -0.29724e-2 ); // a5

        a_coef.push_back(  0.85534e-2 ); // a6
        a_coef.push_back(  0.35322e-2 ); // a7
        a_coef.push_back( -0.48258e-3 ); // a8
        a_coef.push_back( -0.39346e-3 ); // a9
        a_coef.push_back(  0.20338e-3 ); // a10

        a_coef.push_back(  0.25461e-4 ); // a11
        a_coef.push_back( -0.17794e-4 ); // a12
        a_coef.push_back(  0.67394e-5 ); // a13
        a_coef.push_back( -0.21033e-5 ); // a14

        R_cutoff = 8.0;

        delta_Rrms = 0.01; 

        error = 0.01;
    }

    // Calcium 48
    else if ( ( A == 48 ) && ( Z == 20 ) ) {
        a_coef.push_back(  0.44782e-1 ); // a1
        a_coef.push_back(  0.59523e-1 ); // a2
        a_coef.push_back( -0.74148e-2 ); // a3 
        a_coef.push_back( -0.29466e-1 ); // a4
        a_coef.push_back( -0.28350e-3 ); // a5
                                           
        a_coef.push_back(  0.10829e-1 ); // a6
        a_coef.push_back(  0.30465e-2 ); // a7
        a_coef.push_back( -0.10237e-2 ); // a8
        a_coef.push_back( -0.17830e-3 ); // a9
        a_coef.push_back(  0.55391e-4 ); // a10
                                           
        a_coef.push_back( -0.22644e-4 ); // a11
        a_coef.push_back(  0.82671e-5 ); // a12
        a_coef.push_back( -0.27343e-5 ); // a13
        a_coef.push_back(  0.82461e-6 ); // a14
        a_coef.push_back( -0.22780e-6 ); // a15

        R_cutoff = 8.0;

        delta_Rrms = 0.009; 

        error = 0.01;
    }

    // Tin 116
    else if ( ( A == 116 ) && ( Z == 50 ) ) {


        R_coef.push_back( 0.1 ); // R1
        R_coef.push_back( 0.7 ); // R2
        R_coef.push_back( 1.3 ); // R3
        R_coef.push_back( 1.8 ); // R4
        R_coef.push_back( 2.3 ); // R5
        R_coef.push_back( 3.1 ); // R6
        R_coef.push_back( 3.8 ); // R7
        R_coef.push_back( 4.8 ); // R8
        R_coef.push_back( 5.5 ); // R9
        R_coef.push_back( 6.1 ); // R10
        R_coef.push_back( 7.1 ); // R11
        R_coef.push_back( 8.1 ); // R12

        Q_coef.push_back( 0.005727 ); // Q1
        Q_coef.push_back( 0.009643 ); // Q2
        Q_coef.push_back( 0.038209 ); // Q3
        Q_coef.push_back( 0.009466 ); // Q4
        Q_coef.push_back( 0.096665 ); // Q5
        Q_coef.push_back( 0.097840 ); // Q6
        Q_coef.push_back( 0.269373 ); // Q7
        Q_coef.push_back( 0.396671 ); // Q8
        Q_coef.push_back( 0.026390 ); // Q9
        Q_coef.push_back( 0.048157 ); // Q10
        Q_coef.push_back( 0.001367 ); // Q11
        Q_coef.push_back( 0.000509 ); // Q12

        gamma = std::sqrt( 2 / 3.0 ) * 1.60;
        error = 0.03;
        delta_Rrms = 0.001;
        R_cutoff = 12.0;

        for ( unsigned int i = 0; i < R_coef.size(); ++i ) {

            double An = Z * Q_coef[i] 
                          / ( 2 * std::pow( M_PI, 1.5 ) * std::pow( gamma, 3 ) 
                            * ( 1 + 2 * std::pow( R_coef[i] / gamma, 2 ) ) );

            A_coef.push_back( An );
        }
        

    }

    // Tin 124
    else if ( ( A == 124 ) && ( Z == 50 ) ) {

        R_coef.push_back( 0.1 ); // R1
        R_coef.push_back( 0.7 ); // R2
        R_coef.push_back( 1.3 ); // R3
        R_coef.push_back( 1.8 ); // R4
        R_coef.push_back( 2.3 ); // R5
        R_coef.push_back( 3.1 ); // R6
        R_coef.push_back( 3.8 ); // R7
        R_coef.push_back( 4.8 ); // R8
        R_coef.push_back( 5.5 ); // R9
        R_coef.push_back( 6.1 ); // R10
        R_coef.push_back( 7.1 ); // R11
        R_coef.push_back( 8.1 ); // R12

        Q_coef.push_back( 0.004877 ); // Q1
        Q_coef.push_back( 0.010685 ); // Q2
        Q_coef.push_back( 0.030309 ); // Q3
        Q_coef.push_back( 0.015857 ); // Q4
        Q_coef.push_back( 0.088927 ); // Q5
        Q_coef.push_back( 0.091917 ); // Q6
        Q_coef.push_back( 0.257379 ); // Q7
        Q_coef.push_back( 0.401877 ); // Q8
        Q_coef.push_back( 0.053646 ); // Q9
        Q_coef.push_back( 0.043193 ); // Q10
        Q_coef.push_back( 0.001319 ); // Q11
        Q_coef.push_back( 0.000036 ); // Q12

        gamma = std::sqrt( 2 / 3.0 ) * 1.60;
        error = 0.03;
        delta_Rrms = 0.001;
        R_cutoff = 12.0;

        for ( unsigned int i = 0; i < R_coef.size(); ++i ) {

            double An = Z * Q_coef[i] 
                          / ( 2 * std::pow( M_PI, 1.5 ) * std::pow( gamma, 3 ) 
                            * ( 1 + 2 * std::pow( R_coef[i] / gamma, 2 ) ) );

            A_coef.push_back( An );
        }
    }


    // Lead 208
    else if ( ( A == 208 ) && ( Z == 82 ) ) {
        a_coef.push_back(  0.51936e-1 ); // a1
        a_coef.push_back(  0.50768e-1 ); // a2
        a_coef.push_back( -0.39646e-1 ); // a3 
        a_coef.push_back( -0.28218e-1 ); // a4
        a_coef.push_back(  0.28916e-1 ); // a5
                                
        a_coef.push_back(  0.98910e-2 ); // a6
        a_coef.push_back( -0.14388e-1 ); // a7
        a_coef.push_back( -0.98262e-3 ); // a8
        a_coef.push_back(  0.72578e-2 ); // a9
        a_coef.push_back(  0.82318e-3 ); // a10
                                
        a_coef.push_back( -0.14823e-2 ); // a11
        a_coef.push_back(  0.13245e-3 ); // a12
        a_coef.push_back( -0.84345e-4 ); // a13
        a_coef.push_back(  0.48417e-4 ); // a14
        a_coef.push_back( -0.26562e-4 ); // a15
        a_coef.push_back(  0.14035e-4 ); // a16
        a_coef.push_back( -0.71863e-5 ); // a17

        R_cutoff = 12.0;
        delta_Rrms = 0.002; 
        error = 0.01;
    }
    else {

        R_cutoff = 1000;
        std::cout << nucleus_string << " not found." << std::endl;
        return 0;
    }

    
    // r-grid
    double rdelt = 10.0 / 150;
    double rmin = 0.01;
    int rpts = static_cast<int>( R_cutoff / rdelt ) + 1;

    //output file
    std::string filename = "exp_chd_" + nucleus_string 
                         + "_" + param_string + ".dat";
    std::string filename2 = "exp_chd_" + nucleus_string 
                          + "_" + param_string + ".out";

    double norm = 0; // Should be normalized to Ze
    double chd_r2 = 0;
    std::ofstream file( filename.c_str() );
    std::ofstream file2( filename2.c_str() );

    file << rpts << std::endl;
    for( int i = 0; i < rpts; ++i ) {

        double r = rmin + i * rdelt;

        double asum = 0;

        if ( param_string.compare("sog") == 0 ) { // Sum Over Gaussians

            for ( unsigned int n = 0; n < Q_coef.size(); ++n ) {

                asum += A_coef[n]
                    * ( std::exp(-std::pow( ( r - R_coef[n] ) / gamma, 2 ) )
                      + std::exp(-std::pow( ( r + R_coef[n] ) / gamma, 2 ) ) );
            }

        }
        else { // Sum over Spherical Bessels (l = 0)
            for( unsigned int n = 0; n < a_coef.size(); ++n ) {

                asum += a_coef[n] * 
                  boost::math::sph_bessel( 0, ( n + 1 ) * M_PI * r / R_cutoff );
            }
        }

        norm += 4 * M_PI * asum * r * r * rdelt;
        chd_r2 += 4 * M_PI * asum * std::pow( r, 4 ) * rdelt;

        file << r << " " << asum << " " << error * asum << " " 
             << norm << std::endl;

        file2 << r << " " << asum << " " << error * asum << " "
              << norm << std::endl;

    }

    double rms = std::sqrt( chd_r2 / norm );

    //file << rms << " " << delta_Rrms << std::endl;
    std::cout << "norm for " << nucleus_string  << " = " << norm << std::endl;
    std::cout << "RMS radius = " << rms << std::endl;

    file.close();
    file.clear();


    return 1;
}
