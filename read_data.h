#ifndef _READ_DATA_H_
#define _READ_DATA_H_

#include "boost/tuple/tuple.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "store_data.h"
#include "io.h"

data_sets_t
read_volume_integral_data( const std::string &data_dir, 
                           const std::string &prefix,
                           const std::string &suffix ); 
Data_set
read_general_set( const std::string filename ); 

Data_set2
read_matrix ( const std::string filename ); 

Data_set2
read_matrix_ratio ( const std::string filename1, const std::string filename2 ); 

Data_vec 
read_all_data( std::string reactions_filename );

#endif // _READ_DATA_H_
