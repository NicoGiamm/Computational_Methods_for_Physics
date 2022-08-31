#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

extern double pi;

double f (vector<double>& x, int& j, double& k);													//FORZA CH AGISCE SUL PUNTO
double t (double& dt, int& i);														//RESITUISCE IL TEMPO I-ESIMO
double velocity_verlet_X (double& N, double& dt, double& m, double& k, vector<double>& x, vector<double>& v, int& j);	//RESTITUISCE POSIZIONE E VELOCITA' DEL PUNTO DEL SISTEMA
double pot (double& k, double& X);													//ENERGIA POTENZIALE EL SISTEMA
double H (double& X, double& V, double& k, double& m);								//ENERGIA TOTALE DEL SISTEMA

double pot (double& k, double& X){
	
	double potenziale = 0.5*k*pow(X,2);
	
	return potenziale;
}
double H (double&X, double& V, double& k, double& m){
	
	double energia = pow(V,2)/(2*m) + pot(k,X);
	
	return energia;
}

double t (double& dt, int& i){
	
	double t_i = dt*i;
	
	return t_i;
}

double Cn (vector<double>& x, double n){
	
	double Ci = 0;
	 
	for (int j = 1; j < x.size()-1; j++){
		
		Ci += x[j]*sin( (pi*n*j)/(x.size()-1) );
	}
	
	return Ci;
}

double dCn (vector<double>& v, double n){
	
	double dCi = 0;
	
	for (int j = 1; j < v.size()-1; j++){
		
		dCi += v[j]*sin( (pi*n*j)/(v.size()-1) );
	}
	
	return dCi;	
}

double En (double& C, double& dC, int& n, double& k, double& m, int dim){
	
	double Ei = 0;
	
	Ei = (0.5*m)*pow( dC , 2 ) + 2*k*pow( C , 2 )*pow( sin( (pi*n)/(2*(dim+1)) ) , 2 );
	
	return Ei;
}

double f (vector<double>& x, int& j, double& k){
	
	double A = 3;
	double B;
	double forza = k*(x[j+1] - 2*x[j] + x[j-1]) + A*( pow( (x[j+1] - x[j]) , 2 ) - pow( (x[j] - x[j-1]) , 2 ) );
	
	return forza;
}

double velocity_verlet_X (double& N, double& dt, double& m, double& k, vector<double>& x, vector<double>& v, int& j){
	
	double xxx;
	
	xxx = x[j] + v[j]*dt + (1/(2*m))*f(x,j,k)*pow(dt,2);
	
	return xxx;	
}

double velocity_verlet_V (double& N, double& dt, double& m, double& k, vector<double>& x, vector<double>& v, vector<double>& X, int& j){
	
	double vvv;
	
	vvv = v[j] + (1/(2*m))*( f(x,j,k) + f(X,j,k) )*dt;
	
	return vvv;	
}
