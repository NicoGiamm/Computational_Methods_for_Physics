#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

int main(){
	
    double Lx = 150, Ly = 150, x0 = 50, y0 = 75, q = 2, sigma = 12;
    double a = 75, b = 81, V0 = -1.7, Vx;
    double dte = 0.0001;
    double dtl = 0.02;
    int Nx = 310, Ny= 310, Neulero = 200, Nleapfrog = 7500;
    complex<double> pic (0.,1.0);
    
    double hh = Ly/(Ny-1);
    double h = Lx/(Nx-1);
    double norm;
    
    vector<complex<double>*> psi0 (Ny);
    vector<complex<double>*> psiprovv (Ny);
    vector<complex<double>*> psi1 (Ny);
    vector<complex<double>*> psi2 (Ny);
    
 //   ofstream provafile ("lfprova.txt");
	ofstream filecompleto ("leapgenerale.txt");
    
    norm=0.;
    
    for ( int j = 0; j < Ny; j++ ){										
    	
    	psi0[j] = new complex<double>[Nx];
    	psiprovv[j] = new complex<double>[Nx];
    	psi1[j] = new complex<double>[Nx];
    	psi2[j] = new complex<double>[Nx];
    	
    	double y = hh*j;
	
    	for ( int i=0; i < Nx; i++ ){										 
    		
        	double x=h*i;
        
        	psi0[j][i] = psiprovv[j][i] = exp(pic*q*x)*exp( (-pow(x-x0,2)-pow(y-y0,2))/(2.*pow(sigma,2)) );
      		norm+=pow(abs(psi0[j][i]),2);	    
    	}	
    }
     
	norm = norm*Lx*Ly/(Nx*Ny);  
    norm = sqrt(norm);
//    cout << "Norma" <<norm <<'\n';
    
    for( int j = 0; j < Ny; j++ ){
	
    	for(int i = 0; i < Nx; i++){
    		
    		double x = h*i;
    		double y = hh*j;
    		
       		psi0[j][i] = psiprovv[j][i] = psi0[j][i]/norm;
//       		provafile << x << " " << y << " " << pow(abs(psi0[j][i]),2) << endl;
    	}
    }

//CALCOLA PSI1 USANDO EULERO CON UN DT PIU' PICCOLO DI QUELLO CHE SI USERA' IN LEAP FROG MA COMPATIBILE CON LA PRECSIONE DEL METODO FUTURO
    
    for( int n = 0; n < Neulero; n++ ){
    	
    	for ( int j = 1; j < Ny-1; j++ ){
    		
    		for ( int i = 1; i < Nx-1; i++ ){
    			
    			double x = h*i;
    			double y = hh*j;
    			
    			if ( x >= a && x <= b ) Vx = V0;
    			else Vx = 0.;
    			
    			psi1[j][i] = psiprovv[j][i] + (pic*dte/2.)*( ( psiprovv[j][i+1] - 2.*psiprovv[j][i] + psiprovv[j][i-1])/(h*h) +
														( psiprovv[j+1][i] - 2.*psiprovv[j][i] + psiprovv[j-1][i])/(hh*hh) ) - pic*dte*Vx*psiprovv[j][i];
			}
		}
		
		norm = 0.;
		
		for ( int j = 0; j < Ny; j++ ){
    		
    		for ( int i = 0; i < Nx; i++ ){
    			
    			psiprovv[j][i] = psi1[j][i];		
    			norm += pow(abs(psi1[j][i]),2);
    		}
    	}
    	
    	norm = norm*Lx*Ly/(Nx*Ny);
		cout << "Passo " << n << " Norma " << norm << endl;	
	}
	
	//ORA HO PSI0 = T-1, PSI1 = T E POSSO CALCOLARE PSI 2 TRAMITE LEAP-FROG
	//INZIO METODO LEAP-FROG
	
	double varif = 0.;

	for (int n = 0; n < Nleapfrog; n++ ){
		
		for ( int j = 1; j < Ny-1; j++ ){						
    		
    		for ( int i = 1; i < Nx-1; i++ ){
    			
    			double x = h*i;
    			double y = hh*j;
    			
    			if ( x >= a && x <= b ) Vx = V0;
    			else Vx = 0.;
    			
    			psi2[j][i] = psi0[j][i] + (pic*dtl)*( ( psi1[j][i+1] - 2.*psi1[j][i] + psi1[j][i-1])/(h*h) +
														 ( psi1[j+1][i] - 2.*psi1[j][i] + psi1[j-1][i])/(hh*hh) ) - 2.*pic*dtl*Vx*psi1[j][i];
			}
		}
		
		norm = 0.;
		
	//LOOP PER SCRIVERE IL FILE, AGGIORNARE I VETTORI E CALCOLARE LA NORMA AGGIORNATA
		
		for ( int j = 0; j < Ny; j++ ){							 
    		
    		for ( int i = 0; i < Nx; i++ ){
    			
    			double x = h*i;
    			double y = hh*j;
    			
    			psi0[j][i] = psi1[j][i];						//AGGIORNA I VETTORI PER IL CICLO SUCCESSIVO
    			psi1[j][i] = psi2[j][i];
    			norm += pow(abs(psi1[j][i]),2);					//CALCOLA LA NORMA
    			
    			if ( varif == n ){
	        		
	        		filecompleto << x << "\t" << y << "\t" << pow(abs(psi0[j][i]),2) << endl; 		//SCRIVE NEL FILE AD INTERVALLI PREFISATTI
				}
			}
		}
		
		if ( varif == n ){
			filecompleto << endl; 
			varif += 25 ;
		}
		
		norm = norm*Lx*Ly/(Nx*Ny);
		cout << "Passo " << n << " Norma " << norm << endl;		
	}
	
	return 0;		
}
