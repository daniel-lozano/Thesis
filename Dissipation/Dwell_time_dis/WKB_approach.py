import numpy as np
import matplotlib.pyplot as plt



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def k(x,g,Vo,m,ep):#sheck 
    return g*np.sqrt(2*ep*Vo(E-1)/m)-((g*x)**2)/(2*m)


def constant(a,b,c,x,m,hbar):
    return -2*((np.sqrt(2m)/(8*(a**1.5)*hbar))*( 2*np.sqrt(a)*(2*a*x+b)*np.sqrt((a*x+b)*x +c )-(b**2-4*a*c)*np.log(2*np.sqrt(a))*np.sqrt((a*x+b)*x +c )+ 2*a*x+b  )

               
def function(a,b,c,x,m,hbar):
     return -2*((np.sqrt(2m)/(8*(a**1.5)*hbar))*( 2*np.sqrt(a)*(2*a*x+b)*np.sqrt((a*x+b)*x +c )-(b**2-4*a*c)*np.log(2*np.sqrt(a))*np.sqrt((a*x+b)*x +c )+ 2*a*x+b)
                

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio

a=0.5*(20.8)E5#fermi
Nx=1000
x=np.linspace(-a,a,Nx)
dx=abs(x[0]-x[1])

#definiendo energias adimensionales
Ne=100
E=np.linspace(0,2,Ne)

gamma=[0,1E-3,1E-5]# eV/Fermi

#Constantes usadas

Vo=1.8E-0.6 #MeV
m=0.5E6# MeV
hbar=179.3# MeV*fermi

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













