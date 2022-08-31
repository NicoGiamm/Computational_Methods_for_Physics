#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

int main(){
	
    double Lx = 150;
    double dte = 0.0001;
    double dtl = 0.02;
    int Nx = 310, Neulero = 200, Nleapfrog = 50000 ;
    complex<double> pic (0.,1.0);
    
    double h = Lx/(Nx-1);
    double norm;
    
    complex<double>* psi0 = new complex<double>[Nx];
    complex<double>* psiprovv = new complex<double>[Nx];
    complex<double>* psi1 = new complex<double>[Nx];
    complex<double>* psi2 = new complex<double>[Nx];
    
 //   ofstream provafile ("lfprova.txt");
	ofstream filecompleto ("dati_esercizi.txt");
    
    norm=0.;
    
    for ( int j = 0; j < Nx; j++ ){										
 									    		
        double x=h*j;
        if ( x < Lx/2 ){
			psi0[j] = (2/sqrt(Lx))*sin(3.1412*x/Lx);
		}
		else { psi0[j]=0; };
		
      	norm+=pow(abs(psi0[j]),2);	    	
    }
     
	norm = norm*Lx/(Nx);  
    norm = sqrt(norm);
//    cout << "Norma" <<norm <<'\n';
    
    	for(int i = 0; i < Nx; i++){
    		
    		double x = h*i;
    		
       		psi0[i] = psiprovv[i] = psi0[i]/norm;
//       		provafile << x << " " << y << " " << pow(abs(psi0[j][i]),2) << endl;
    	}
    

//CALCOLA PSI1 USANDO EULERO CON UN DT PIU' PICCOLO DI QUELLO CHE SI USERA' IN LEAP FROG MA COMPATIBILE CON LA PRECSIONE DEL METODO FUTURO
    
    for( int n = 0; n < Neulero; n++ ){
    		
		for ( int i = 1; i < Nx-1; i++ ){
    			
    		double x = h*i;	
    		psi1[i] = psiprovv[i] + (pic*dte/2.)*( psiprovv[i+1] - 2.*psiprovv[i] + psiprovv[i-1])/(h*h);								
		}
		
		norm = 0.;
		
  		for ( int i = 0; i < Nx; i++ ){  			
    			psiprovv[i] = psi1[i];		
    			norm += pow(abs(psi1[i]),2);
    		}
    	
    	norm = norm*Lx/(Nx);
		cout << "Passo " << n << " Norma " << norm << endl;	
	}
	
	//ORA HO PSI0 = T-1, PSI1 = T E POSSO CALCOLARE PSI 2 TRAMITE LEAP-FROG
	//INZIO METODO LEAP-FROG
	
	double varif = 0.;

	for (int n = 0; n < Nleapfrog; n++ ){
					
  		for ( int i = 1; i < Nx-1; i++ ){  	
					
    		double x = h*i;			
    		psi2[i] = psi0[i] + (pic*dtl)* ( psi1[i+1] - 2.*psi1[i] + psi1[i-1])/(h*h) ;
    		
		}
														 
		norm = 0.;
		
	//LOOP PER SCRIVERE IL FILE, AGGIORNARE I VETTORI E CALCOLARE LA NORMA AGGIORNATA
								 	
    		for ( int i = 0; i < Nx; i++ ){
    			
    			double x = h*i;
    			
    			psi0[i] = psi1[i];						//AGGIORNA I VETTORI PER IL CICLO SUCCESSIVO
    			psi1[i] = psi2[i];
    			norm += pow(abs(psi1[i]),2);					//CALCOLA LA NORMA
    			
    			if ( varif == n ){
	        		
	        		filecompleto << x << "\t" << pow(abs(psi0[i]),2) << endl; 		//SCRIVE NEL FILE AD INTERVALLI PREFISATTI
				}
			}
		
		if ( varif == n ){
			filecompleto << endl; 
			varif += 50 ;
		}
		
		norm = norm*Lx/Nx;
		cout << "Passo " << n << " Norma " << norm << endl;		
	}
	
	return 0;		
}
