#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "velocityverlet.h"
#include "Energia.h"

using namespace std;

double pi = 3.14159;

double t (double& dt, int& i);

int main(){
	
	int N = 86000;
	int num = 33;	
	double m = 1;
	double k = 1;
	double A = 0.1;	
	double dt = 0.1;
	int n = 1;
		
	Campo campo(m,k,A,dt,num);
	Energia energia(m,k);
	vector<double> E(5,0);

	ofstream file1 ("E1.txt");
	ofstream file2 ("E2.txt");
	ofstream file3 ("E3.txt");
	ofstream file4 ("E4.txt");
	ofstream file5 ("E5.txt");
		
	for (int i = 0; i < N ; i++){
		
		campo.evolve();
		
		for (n = 1; n <= E.size(); n++){
			
			E[n-1] = energia.compute( n , campo );
		}
		
		file1 << t(dt,i) << "\t" << E[0] << endl;
		file2 << t(dt,i) << "\t" << E[1] << endl;
		file3 << t(dt,i) << "\t" << E[2] << endl;
		file4 << t(dt,i) << "\t" << E[3] << endl;
		file5 << t(dt,i) << "\t" << E[4] << endl;
		
		campo.x = campo.X;
		campo.v = campo.V;
	}
	
	cout << "fatto" << endl;
	
	return 0;
}

double t (double& dt, int& i){
	
	double t_i = dt*i;
	
	return t_i;
}											
