/* String manipulation and I/O utilities.
 * 
 * Mark Burnett, November 2008
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "io.h"

namespace util {

// --------------------------------------------------------------------
// String related utilities
// --------------------------------------------------------------------
std::vector< std::string >
split(const std::string &s) {
    std::vector< std::string > v;
    std::string my_copy(s);
    boost::algorithm::trim(my_copy);
    boost::split(v, my_copy, boost::is_any_of(" \t"),
            boost::algorithm::token_compress_on);
    return v;
}

std::string IntToStr( int n ) {

    std::ostringstream result;
    result << n;
    return result.str();

}

// --------------------------------------------------------------------
// I/O related utilities
// --------------------------------------------------------------------
// Reads lines from commented files and puts them in a list.
std::list< std::string >
read_commented_file(const std::string &filename) {
    std::ifstream file( filename.c_str() );
    if (!file.is_open()) {
        std::cout << "file " << filename << " (could not be opened) " << std::endl; 
        throw std::exception();
    }

    std::string line;
    std::list< std::string > result;
    while (std::getline(file, line) ) {
        if ( !line.empty() && '#' != line[0] )
        result.push_back(line);
    }

    return result;
}

} // end namespace util
