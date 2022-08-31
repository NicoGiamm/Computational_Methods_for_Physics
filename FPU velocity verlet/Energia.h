#include <iostream>
#include <cmath>
#include <vector>
#include "velocityverlet.h"

#ifndef Energia_h
#define Energia_h

using namespace std;

class Energia {
	
	public:
		
		Energia( double& m , double& k );
		~Energia();
			
		double compute( int& n , Campo& Cv );
		
	private:
		
		double C_;
		double dC_;
		double E_;
		double m_, k_;
		
		void Cn( int& n , Campo& Cv );
		void dCn( int& n , Campo& Cv );
		void En( int& n , Campo& Cv );
};

#endif
