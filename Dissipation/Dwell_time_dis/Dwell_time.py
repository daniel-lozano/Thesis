import numpy as np
import matplotlib.pyplot as plt



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def Calc_k(m,Vo,hbar,E): # k antes de entrar a la barrera

    return np.sqrt(2*m*Vo*E)/hbar


def Calc_kp(m,Vo,hbar,ep,E,gamma,x):# k  dentro de la barrera
    
    return np.sqrt(2*m*Vo*ep*(E-1))/hbar

    
def Calc_kt(m,Vo,hbar,ep,E,gamma,L):# k  afuera de la barrera asumiendo una perdida repentina de energia

    termino1=gamma*L*np.sqrt(2*ep*(E-1)/(m*Vo))
    
    termino2=(ep*(gamma*L)**2)/(2*m*Vo)
    
    return (np.sqrt(2*m*Vo)/hbar)*np.sqrt(ep*(E-1-termino1+termino2))


def Calc_F2(k,kp,kt,a): #calculo del coeficiente de transmision

    return ( (2*kp*kt)**2 )/( ((kp*k+kt*kp)*np.cos(2*kp*a))**2  + ( (kp**2+kt*k)*np.sin(2*kp*a))**2 )

def Calc_C2(F,k,kp,kt): #calculo del coeficiente C dentro de la barrera

    return 0.25*F*(1+kt/kp)**2

def Calc_D2(F,k,kp,kt): #calculo del coeficiente D dentro de la barrera
    
    return 0.25*F*(1-kt/kp)**2

def Calc_CD_DC(F,k,kp,kt,a):#calculo de C D^{*}

    return 0.5*F*(1-kt/kp)*(1+kt/kp)*np.cos(2*kp*a)

def p(m,Vo,hbar,ep,E,gamma,x):# k  afuera de la barrera asumiendo una perdida repentina de energia

    termino1=gamma*x*np.sqrt(2*ep*(E-1)/(m*Vo))
    
    termino2=(ep*(gamma*x)**2)/(2*m*Vo)
    
    return (np.sqrt(2*m*Vo)/hbar)*np.sqrt(ep*(E-1-termino1+termino2))

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio

a=10 #Fermi
Nx=1000
x=np.linspace(-a,a,Nx)
dx=abs(x[0]-x[1])

#definiendo energias adimensionales
Ne=100
E=np.linspace(0,2,Ne)

gamma=[0,1E-3,1E-5]# eV/Fermi

#Constantes usadas

Vo=1.8 #eV
m=0.5E2#E9 #eV
hbar=6.582119514E-16 # eV*second

#matriz que traza el cambio entre E/Vo
epsilon=np.zeros(len(E))
#aqui se guardaran los tiempos relacionados con la energia
T=np.zeros(len(E))

for i in range(Ne):
    if(E[i]<1):
        epsilon[i]=-1
    if(E[i]>1):
        epsilon[i]=1

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comenzando el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


for i in range(len(E)):

    suma=0
    
    k=Calc_k(m,Vo,hbar,E[i])
    kp=Calc_kp(m,Vo,hbar,epsilon[i],E[i],gamma[0],x)
    kt=Calc_kt(m,Vo,hbar,epsilon[i],E[i],gamma[0],2*a)

    P=p(m,Vo,hbar,epsilon[i],E[i],gamma[0],x)
    
    F=Calc_F2(k,kp,kt,a)
    #print("transmission",F,"index",i)#hay algo raro aqui
    
    C=Calc_C2(F,k,kp,kt)
    D=Calc_D2(F,k,kp,kt)
    
    CD_DC=Calc_CD_DC(F,k,kp,kt,a)

    
    for j in range(len(x)):
        '''
        termino1=C*np.exp(2*P[j]*x[j])
        
        termino2=D*np.exp(-2*P[j]*x[j])
        
        termino3=CD_DC
        
        suma+=dx*(termino1+termino2+termino3)
        '''
        suma+=i
    T[i]=suma/(hbar*k*F/m)













