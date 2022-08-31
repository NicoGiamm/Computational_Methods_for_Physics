#include <iostream>
#include <cmath>
#include <vector>
#include "Energia.h"

using namespace std;

extern double pi;

Energia::Energia( double& m , double& k ):m_(m), k_(k){
}

Energia::~Energia(){
}

void Energia::Cn( int& n, Campo& Cv ){
	
	C_ = 0;
	
	for (int j = 1; j < Cv.x.size()-1; j++){
		
		C_ += Cv.x[j]*sin( (pi*n*j)/(Cv.x.size()-1) );
	}
	
	return;
}

void Energia::dCn( int& n , Campo& Cv ){
	
	dC_ = 0;
	
	for (int j = 1; j < Cv.v.size()-1; j++){
		
		dC_ += Cv.v[j]*sin( (pi*n*j)/(Cv.v.size()-1) );
	}
	
	return;
}

void Energia::En( int& n , Campo& Cv ){
	
	E_ = 0;
	
	E_ = (0.5*m_)*pow( dC_ , 2 ) + 2*k_*pow( C_ , 2 )*pow( sin( (pi*n)/(2*(Cv.x.size()-1)) ) , 2 );
	
	return;
}

double Energia::compute( int& n , Campo& Cv ){
	
	Energia::Cn( n , Cv );
	Energia::dCn( n , Cv );
	Energia::En( n , Cv );
	
	return E_;
}
