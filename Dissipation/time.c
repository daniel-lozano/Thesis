#include <stdio.h>    
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "inicial.h"
#define e 2.71828182845904523536
#define PI 3.14159265358979323846

double integrar(double dx,double *x,double *fx,int X);

void dissfunc(double *gamma,double *t,double *Ft,int Time);

void iniciarOMG(double *omega);

void iniciarT(double *t);

void t_func_dis_above(double complex *w,double E,double Vop, double complex *Vo,double d, double m, double hbarc,int SIZE,int OMG);

void t_func_dis_bellow(double complex *w,double E,double Vop, double complex *Vo,double d, double m, double hbarc,int SIZE,int OMG);

double complex fourier(double time, double complex *w, double *freq, double df, int OMG);

double integrate_n_square(double dt, double complex *w, double *ft, double TIME,int shift);


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(){

  
 
  int i,j,k;
  int Time=1000;
  int OMG=10000;
  double suma=0;
  double dt=1.0/Time;
  double dw=1.0/OMG;
  double *omega;
  double complex *Vo; //frecuencias posibles y arreglo de potencial 
  double *t; //tiempo usado
  double *Ft;//funciones de disipacion
  double *gamma; //gammas posibles

  omega=malloc(sizeof(double)*OMG);//arreglo de omega 

  Vo=malloc(sizeof(double)*OMG);

  t=malloc(sizeof(double)*Time);//arreglo de tiempo 

  Ft=malloc(sizeof(double)*Time*2);//arreglo de disipacion

  gamma=malloc(sizeof(double)*2);//arreglo de gammas



  //iniciando arreglos de frecuencias y tiempo--------------------------------------------

  iniciarOMG(omega);
  iniciarT(t);
  
  printf("valor min omega= %f valor max omega= %f\n",omega[0],omega[OMG]);
  printf("valor min t= %f valor max t= %f\n",t[0],t[Time]);

  
  //Definiendo gammas posibles------------------------------------------------------------
  
  gamma[0]=1E-2;
  gamma[1]=5E-2;

  //mirando la funcion de disipacion------------------------------------------------------

  dissfunc(gamma,t,Ft,Time);

  Ft[0]=1.0;
  Ft[Time]=1.0;
  
  for(i=0;i<2*Time;i++){
    Ft[i]=sqrt(Ft[i]);
  }

  FILE *f;

  f=fopen("dissipation_function_sqrt.dat","w");

  for(i=0;i<=Time;i++){
    fprintf(f,"%f %f %f\n",t[i],Ft[i],Ft[i+Time]);
  }
  fclose(f);

 //creando variables y arrelos que se usaran-------------------------------------------

  double Vop=1.8; // eV
  double hbarc=0.19733; // ev Microm *c
  double c=2.998E14;// Microm/sec
  double hbar=0.065821E-14;
  double m=0.511E6; // ev
  double d=20.8E-4; // Microm
  
  printf("\nVop=%f eV\n",Vop);
  printf("\nhbarc=%f eV*microM*c\n",hbarc);
  printf("\nc=%f microM/sec\n",c);
  printf("\nhbar=%f eV*microM\n",hbar);
  printf("\nm=%f eV\n",m);
  printf("\nd=%f microM\n",d);
  printf("prueba %lf \n",hbar);
  double *E;
  int steps=200;
  double de=1.0/steps;
  
  E=malloc(sizeof(double)*steps);
  
  for(int i=0;i<steps;i++){
    E[i]=2*de*(i+1);
    //printf("E[%d]=%f\n",i,E[i]);
  }
  
  for(i=0;i<OMG;i++){
    Vo[i]=Vop-hbar*omega[i];
    //printf("Vo[%d]=%f\n",i,Vo[i]);
  
  }

  
  double complex *w;
  double complex *FuncD;
  FuncD=malloc(sizeof(double complex)*OMG);
  w=malloc(sizeof(double complex)*Time);

  
  
  double Energy=0.9;
  for(j=0;j<Time;j++){
    t_func_dis_above(FuncD,Energy,Vop,Vo,d,m,hbarc,0,OMG);
    
    w[j]=Ft[j]*fourier(t[j],FuncD,omega,dw,OMG);
  
  }


  

  return 0;
}

//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------
//--------------------------------------------------------

void iniciarT(double *t){
  int i;
  int Time=1000;
  double dt=0.01;
  for(i=0;i<=Time;i++){
    t[i]=i*dt;
  }
}

void iniciarOMG(double *omega){
  int i;
  int OMG=10000;
  double dw=0.1;
  for(i=0;i<=OMG;i++){
    omega[i]=-1000+2.0*i*dw;
  }
}

double integrar(double dx,double *x,double *fx,int X){
  double suma=0.0;
  int i;

  for(i=0;i<X;i++){
    suma+=dx*fx[i];
    }
  return suma;
}

 void dissfunc(double *gamma,double *t,double *Ft,int Time){

  double OMEGA= 100.0;
  double sigma=0;
  double g=0.0;
  int i,j;
  for(j=0;j<2;j++){
    g=gamma[j];
    sigma=sqrt(g*g+pow(4.0*g*OMEGA/PI,2));
  
    for(i=0;i<=Time;i++){
      Ft[i+Time*j]=sigma*t[i]*exp(g*t[i])/sinh(sigma*t[i]);
    }
  }

}


void t_func_dis_above(double complex *w,double E,double Vop, double complex *Vo,double d, double m, double hbarc,int SIZE,int OMG){
  
  int i;
  double complex im=I*1.0;
  double k;
  double complex kp;

  for(int i=0;i<OMG;i++){
    k=sqrt(2*m*E*Vop)/hbarc;
    kp=sqrt(2*m*(E-1.0)*Vo[i])/hbarc;
  
    w[i]=-(2.0*im*k*kp*cexp(-im*k*d)) / ( (k*k + kp*kp)*csin(kp*d) + 2*im*k*kp*ccos(kp*d) ) ;
   
  }

}

void t_func_dis_bellow(double complex *w,double E,double Vop, double complex *Vo,double d, double m, double hbarc,int SIZE,int OMG){
  
  int i;
  double complex im=I*1.0;
  double k;
  double complex kp;

  for(int i=0;i<OMG;i++){
    k=sqrt(2*m*E*Vop)/hbarc;
    kp=sqrt(2*m*(1.0-E)*Vo[i])/hbarc;
  
    w[i]=(2.0*k*kp*cexp(-im*k*d)) / ( (k*k - kp*kp)*im*csinh(kp*d) - 2*k*kp*ccosh(kp*d) );
    
   
  }

}

double integrate_n_square(double dt, double complex *w, double *ft, double TIME,int shift){
  int i;
  double complex suma=0;
  double total=0;

  for(i=0;i<TIME;i++){
    suma+=dt*w[i]*ft[i+shift];
  }
  total=creal(suma)*creal(suma)+cimag(suma)*cimag(suma);

  return total;
}

double complex fourier(double time, double complex *w, double *freq, double df, int OMG){
  
  int i;
  double complex suma=0;

  for(int i=0; i<OMG;i++){
    suma+=df*cexp(-I*time*freq[i])*w[i]/(2*PI);
  }
  
  return suma;
}
