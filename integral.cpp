#include<stdio.h>
#include<iostream>
#include<cmath>
#include "boundRspace.h"

main(){

 std::ifstream file("/scratch/mmahzoon/nonlocal/r_space_scattering/waves/Results_DOM/Spectral/spectral_p_u.out");

 double energy;
 double spectral; 

 std::vector<double> a;
 std::vector<double> b;
 int line = 327;
 for (int ii=0; ii < line ; ++ii){
	file >> energy >> spectral; 
        a.push_back(energy);
        b.push_back(spectral);	
 }
 double sum = 0.0;
 int jj=0;

 while( a[jj] < -0.0){
	sum = sum + b[jj] * (a[jj+1]-a[jj]);	
	jj = jj+1; 
 }
file.close();
cout << sum <<std::endl;
return 0;
}

