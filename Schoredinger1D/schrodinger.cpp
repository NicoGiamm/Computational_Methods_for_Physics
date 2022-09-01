#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>

using namespace std;

void solve_tridiagonal(int n, complex<double> *d, complex<double> *u, complex<double> *l,complex<double> *b, complex<double> *a ){
  //soluzione del problema Ma=b con M matrice tridiagonale e b dato
  //n-dimensione matrice, d parte diagonale, u diagonale superiore,l diagonale inferiore, vettore b, a soluzione (output)


    
    complex<double>* alfa = new complex<double>[n];
    complex<double>* beta = new complex<double>[n];
    
    alfa[0]=-d[0]/u[0];
    beta[0]=b[0]/u[0];
    
    for(int i=1;i < n; i++){
        alfa[i]=(-l[i-1]/(u[i]*alfa[i-1])-d[i]/u[i]);
        beta[i]=b[i]/u[i]+l[i-1]*beta[i-1]/(u[i]*alfa[i-1]);
                 
    }
    
    
    a[n-1]=(b[n-1]/l[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+d[n-1]/l[n-2]);
    
    for(int i=n-1;i>0;i--){
        a[i-1]=a[i]/alfa[i-1]-beta[i-1]/alfa[i-1];
    }
    
    
    free(alfa);
    free(beta);
    
    
}


int main (){

	double planck = 1;
	double m = 1;
	double pi = acos(-1);
	double L = 500;
	double x0 = 200;
	double q = 2;
	double sigma = 20;
	double a = 250;
	double b = 260;
	double V0 = 1.7;
	double Vi;
	int Nx = 2000;
	double Nt = 2000;
	double dt = 0.1;
	double dx = L/Nx;
	double xi;
	complex<double> imm (0.,1.0);
	complex<double>* psi0 = new complex<double>[Nx];
	complex<double>* psi = new complex<double>[Nx];
	complex<double>* d = new complex<double>[Nx];
	complex<double>* u = new complex<double>[Nx];
	complex<double>* l = new complex<double>[Nx];
	
	ofstream file ("schrodingaaa.txt");

	for (int i = 0; i < Nx; i++){
		
		xi = i*dx;

		psi0[i] = exp( imm*q*xi ) * exp( -( pow( xi-x0 , 2 )/( 2*pow( sigma , 2 ) ) ) );	
		
		file << xi << "	" << pow( abs(psi0[i]) , 2 ) << endl; 
	}

	u[0] = 1.;
	d[0] = d[Nx-1] = imm*( 4*m*pow(dx,2)/(planck*dt) ) - 2. - ( 2*m*pow( dx , 2 )/(pow( planck , 2 ) ) )*1000000.;
	l[Nx-1] = 1.;
	
	
	for ( int i = 1; i < Nx-1; i++ ){
		
		xi = i*dx;

		if ( xi <= 260 && xi >= 250 ) Vi = V0;
		else Vi = 0;
		
		d[i] = imm*( 4*m*pow(dx,2)/(planck*dt) ) - 2. - ( 2*m*pow( dx , 2 )/(pow( planck , 2 ) ) )*Vi;
		u[i] = 1.;
		l[i] = 1.;
	}
	
	for ( int i = 0; i < Nt; i++ ){

		solve_tridiagonal( Nx , d , u , l , psi0 , psi );

		for ( int j = 0; j < Nx; j++ ){

			xi = j*dx;

			file << xi << "	" << pow( abs(psi[j]) , 2 ) << endl; 

			psi0[j] = psi[j];
		}
	} 
	

	return 0;
}
