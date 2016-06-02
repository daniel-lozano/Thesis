#include <stdio.h>    
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define e 2.71828182845904523536
#define PI 3.14159265358979323846

double integrar(double dx,double *x,double *fx,int X);
void dissfunc(double *gamma,double *t,double *Ft,int Time);
void iniciarOMG(double *omega);
void iniciarT(double *t);


int main(){

 
  int i;
  int Time=1000;
  int OMG=10000;
  double suma=0;
  double complex *x;
  double *omega; //frecuencias posibles
  double *t; //tiempo usado
  double *Ft;//funciones de disipacion
  double *gamma; //gammas posibles

  omega=malloc(sizeof(double)*OMG);//arreglo de omega 
  t=malloc(sizeof(double)*Time);//arreglo de tiempo 
  Ft=malloc(sizeof(double)*Time*2);//arreglo de disipacion
  gamma=malloc(sizeof(double)*2);//arreglo de gammas
  x=malloc(sizeof(double complex)*1);//arreglo complejo


  x[0]=(double complex)I;
  printf("valor de x %f+ i%f\n",creal(x[0]),cimag(x[0]));

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

  FILE *f;

  f=fopen("dissipation.dat","w");

  for(i=0;i<=Time;i++){
    fprintf(f,"%f %f %f\n",t[i],Ft[i],Ft[i+Time]);
  }
  fclose(f);

  
 


  return 0;
}


//---------------Functions----------------//

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
