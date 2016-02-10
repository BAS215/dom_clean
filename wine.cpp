#include <cmath>
#include<iostream>
// #include "boundRspace.h"
#include<vector>
#include<fstream>

int main() {
std::vector <double> r;
std::vector <double> V;

double rn=10;
double r0=.01;
int Npoints=1000;
double dr;
dr=(rn-r0)/Npoints;
double V_ws=-202.21;
double R_ws=0.93;
double beta=0.87;
double a_ws=0.66;
double rho=0.81;

std::ofstream winein;
winein.open("wine_4He_potential.out");

for (int i=0;i<Npoints;++i){
	r.push_back(r0+i * dr);
	V.push_back(V_ws * (1/(1+std::exp((r[i]-R_ws)/a_ws))-beta * std::exp(-std::pow(r[i]/rho,2))));
	winein << r[i] << " " << V[i] << " " << 0 << std::endl;	
    }	


 return 0;

}	
