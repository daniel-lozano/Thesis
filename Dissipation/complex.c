#include <stdio.h>    
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define e 2.71828182845904523536
#define PI 3.14159265358979323846

int main(){
  double complex a0;
  double complex a1;
  double complex a2;
 

  //mirando las funciones de exponenciacion

  a0=cexp(I*2*PI);// como es complejo se debe usar la funcion cexp y no exp
  a1=cexp(I*PI/2.0);

  
  a2=a0+a1;
  printf("numero complejo a2, im(a2)=%f,re(a2)=%f\n",cimag(a2),creal(a2));

  //mirando las funciones de fracciones 

  a0=1/I;
  a1=(1/I)*I;
  a2=a0+a1;
  
  a2=a0+a1;
  printf("numero complejo a2, im(a2)=%f,re(a2)=%f\n",cimag(a2),creal(a2));
  
  double complex *A;
  double complex *B;
  double  *C;
  int N=15;
  int i;
  
  A=malloc(sizeof(double)*2);
  A[0]=1+I;
  A[1]=a1;
  printf("numero complejo A[0], im(A[0])=%f,re(A[0])=%f\n",cimag(A[0]),creal(A[0]));
  
  B=malloc(sizeof(double)*N);
  C=malloc(sizeof(double)*2*N);
  double a,b;
  double complex c;

  for(i=0;i<N;i++){
    srand48(i*i);
    a=(0.5-drand48())*2;
    b=(0.5-drand48())*2;
    c=a+b*I;
    B[i]=c;
    C[2*i]=creal(c);
    C[(2*i+1)]=cimag(c);
    printf("re(B[%d])=%f,im(B[%d])=%f,\n",i,creal(B[i]),i,cimag(B[i]));
    printf("re(C[%d])=%f,im(C[%d])=%f,\n",2*i,C[2*i],2*i+1,C[(2*i+1)]);
    
    
  }
  
  



  return 0;
  

  
}
