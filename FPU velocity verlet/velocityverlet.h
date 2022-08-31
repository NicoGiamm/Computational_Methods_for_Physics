#include <iostream>
#include <cmath>
#include <vector>

#ifndef velocityverlet_h
#define velocityverlet_h

using namespace std;

class Campo {
	
	public:
		
		Campo( double& m, double&k, double& A, double& dt, int& num );
		~Campo();
			
		vector<double> x;
		vector<double> v;
		vector<double> X;
		vector<double> V;
		
		double f( int& j, vector<double> x_);
		void evolve ();
		
	private:
		
		const double dt_;
		const double m_;
		const double k_;
		const double A_;
		const double num_;
};

#endif
