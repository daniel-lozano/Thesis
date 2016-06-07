#include <stdio.h>    
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define e 2.71828182845904523536
#define PI 3.14159265358979323846


int main(){
  double *x;
  double *w;
  double *y;
  double *wy;
  int i,j;
  int N=100;
  int M=1000;
  double dx=1.0/N;
  double dw=0.01;
  FILE *f,*g;

  f=fopen("gaussian.dat","w");
  g=fopen("fourier_gaussian.dat","w");

  x=malloc(sizeof(double)*N);
  y=malloc(sizeof(double)*N);

  w=malloc(sizeof(double)*M);
  wy=malloc(sizeof(double)*M);

  //llenando el archivo de gaussiana
  for(i=0;i<N;i++){
    x[i]=-10+20*i*dx;
    y[i]=exp(-pow(x[i],2));
    fprintf(f,"%f %f\n",x[i],y[i]);
    
    }

  fclose(f);
  
  for(i=0;i<M;i++){
    w[i]=-10+2*i*dw;
  }
  printf("w[0]=%f w[M-1]=%f\n",w[0],w[M-1]);


  double complex suma;  

  for(j=0;j<M;j++){

    suma=0;

    for(i=0;i<N;i++){

      suma+=dx*cexp(-I*x[i]*w[j])*y[i]/(2*PI);
    
    }
    
    wy[j]=creal(suma);//sqrt(pow(creal(suma),2)+pow(cimag(suma),2));
    fprintf(g,"%f %f\n",w[j],wy[j]);

  }
  fclose(g);
  
  
  
  
  
  



  return 0;
}
