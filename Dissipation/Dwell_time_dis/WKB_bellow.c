#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define e 2.71828182845904523536
#define PI 3.14159265358979323846




//%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


double k(double m,double eta,double Vo,double E,double hbar,double x);

double function(double m,double eta,double Vo,double E,double hbar,double x);


int main(){

    int i,j,l=0;
    
    
    //definiendo energias adimensionales
    int Ne=100;
    int neta=5;
    int Nx=10;
    double *E,*eta,*Time,*x;
    
    E=malloc(sizeof(double)*Ne);
    Time=malloc(sizeof(double)*Ne);
    eta=malloc(sizeof(double)*neta);
    x=malloc(sizeof(double)*Nx);

    
    
    //Constantes usadas
    
    double a=20.8E5; //fermi
    double Vo=1.8E-6; // MeV
    double m=0.5E6; // MeV
    double hbar=197.327;// MeV*fermi =[hbar*c]
    double dx=a/Nx;
    
    
    
    printf("\nmass=%f  MeV\n",m);
    printf("L=%f fermi\n",a);
    printf("Vo=%f MeV\n",Vo);
    printf("hbar=%f Mev*fermi\n",hbar);

    //llenando el arreglo de Energias----------------------------------------------
    
    for(i=0;i<Ne;i++){
    
        E[i]=0.005+(1.0/Ne)*i;
        
    }
    printf("\nEnergy in [%f,%f]\n",E[0],E[Ne-1]);
    
    double hbar_k_l=sqrt(2.0*m*Vo*(1-E[Ne-1]))/a;
    printf("\nCondicion:\neta<< %.20f \n",hbar_k_l);
    
    //llenando arreglo de valores de dissipacion----------------------------------------------

    for(i=0;i<neta;i++){
        eta[i]=hbar_k_l * pow(10,-i-2);
        printf("eta[%d]=%.15f\n",i,eta[i]);
    
    }
    
    //llenando arreglo de valores de posicion----------------------------------------------

    for(i=0;i<Nx;i++){
        x[i]=(i)*a*1.0/Nx;
        
    }
    printf("\nx[%d]=%.15f\n",0,x[0]);
    printf("x[%d]=%.15f\n",Nx-1,x[Nx-1]);
    
    
    //mirando los posibles valores de k
    
    for(i=0;i<Ne;i++){
        double suma=0;
        
        for(j=0;j<Nx;j++){
    
            printf("probando %f\n",function(m,eta[1],Vo,E[i],hbar,x[j]));
            
            suma+=(m/hbar)*dx*k(m,eta[1],Vo,E[i],hbar,x[j]);
            }
        Time[i]=suma;
    }
    
    
    
    
    
    
    
    return 0;
}

//funcitones que sera utilizadas

double k(double m,double eta,double Vo,double E,double hbar,double x){
    
    float k0=sqrt(2*m*Vo*(1-E))/hbar;
    
    return sqrt(k0*k0+(2*eta*k0*x)/hbar +pow(eta*x/hbar,2)) ;
}

double function(double m,double eta,double Vo,double E,double hbar,double x){
    
    double k0=sqrt(2*m*Vo*(1-E))/hbar;
    
    double termino1=-0.5*(2*x+2*hbar*k0/eta);
    
    double termino2=sqrt(pow(k0,2) +(2*eta*k0*x)/hbar +pow(eta*x/hbar,2)  ) ;
    
    double termino3= (1/eta)*(4*hbar*pow(k0,2));
    
    return (termino1*termino2+termino3);
}
