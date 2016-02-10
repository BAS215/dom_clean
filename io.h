/* String manipulation and general I/O stuff.
 * 
 * Mark Burnett, November 2008
 */

#ifndef _UTIL_IO_H_
#define _UTIL_IO_H_

#include <vector>
#include <list>
#include <string>
#include <sstream> // SJW 08/09/2010
#include <fstream>
#include <stdexcept>
#include <iostream>

namespace util {

// --------------------------------------------------------------------
// Provide standard split functionality for c++, using boost's weird split.
// --------------------------------------------------------------------
std::vector< std::string > split(const std::string &s);

std::string IntToStr( int n ); //SJW 08/09/2010
// --------------------------------------------------------------------
// I/O related utilities
// --------------------------------------------------------------------
std::list< std::string > read_commented_file(const std::string &filename);

} // end namespace util


#endif // _UTIL_IO_H_
