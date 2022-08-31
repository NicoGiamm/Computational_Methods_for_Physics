#include <iostream>
#include <cmath>
#include <vector>
#include "velocityverlet.h"

using namespace std;

extern double pi;

Campo::Campo( double& m, double&k, double& A, double& dt, int& num ): dt_(dt), m_(m), k_(k), A_(A), num_(num){

	for (int j = 0; j < num_; j++){
	
		x.push_back(0);
		X.push_back(0);
		v.push_back(0);
		V.push_back(0);
	}
	
	for (int j = 1; j < num_-1; j++){
		
		x[j] =  A_*sin( (pi*j)/(num_-1));
	}
}

Campo::~Campo(){
}

double Campo::f( int& j, vector<double> x_ ){
	
	double alpha = 3;
	
	double forza = k_*(x_[j+1] - 2*x_[j] + x_[j-1]) + alpha*( pow( (x_[j+1] - x_[j]) , 2 ) - pow( (x_[j] - x_[j-1]) , 2 ) );
	
	return forza;	
}

void Campo::evolve(){
	
	for ( int j = 1; j < x.size()-1; j++ ){
		
		X[j] = x[j] + v[j]*dt_ + (1/(2*m_))*Campo::f(j,x)*pow(dt_,2);
	}
	
	for ( int j = 1; j < v.size()-1; j++ ){
		
		V[j] = v[j] + (1/(2*m_))*( Campo::f(j,x) + Campo::f(j,X) )*dt_;
	}
	
	return;
}
