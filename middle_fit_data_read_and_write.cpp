//this file read a file ("MiddleFit.inp") of one line of the  fitted data in the middle of the fit and makes the new inputfileas "hosfitMiddleFit.inp" to 
//feed to the "domgen" or "chisq" or anything else.
//(fitted dataline must be rescaled as well)
//
//
//
//For simplicity we take into account the file 
//
//
//
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include "types.h"
#include "io.h"
#include <list>

int main(){
     std::vector<double> unscaled_input_numbers;
     double input;
    //reading Input file that is used to fit the data
    //
    
//    std::ifstream fitted_data_line("hosfitApril131804.inp");
    
    //reading the fitted data file line which is not scaled

    std::ifstream fitted_data_line("MiddleFit.inp");
    while(!fitted_data_line.eof()) 
    { 
    //Get input 
     fitted_data_line >> input;
     if(fitted_data_line.eof() ) break;
     unscaled_input_numbers.push_back(input);
    ////Print input 
    //std::cout << input<<std::endl; 
    } 

    for (unsigned int k=0;k< unscaled_input_numbers.size() ; ++k){
        std::cout<< unscaled_input_numbers[k]<<std::endl;
    }


    //
    ////Close file 
    fitted_data_line.close(); 
    //
    return 0; 
} 
    





