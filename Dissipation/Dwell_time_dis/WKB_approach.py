import numpy as np
import matplotlib.pyplot as plt



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def k(m,eta,Vo,E,hbar,ep):#sheck
    
    k0=np.sqrt(2*m*ep*Vo*(E-1))/hbar
    
    return np.sqrt(k0**2-(ep*2*eta*k0*x)/hbar +(eta*x/hbar)**2 )

def function(m,eta,Vo,E,hbar,ep):
    
    k0=np.sqrt(2*m*ep*Vo*(E-1))/hbar
    
    termino1=-0.5*(2*x-ep*hbar*k0/eta)
    
    termino2=np.sqrt(k0**2 -(ep*2*eta*k0*x)/hbar +(eta*x/hbar)**2   )
    
    termino3= (1/eta)*(np.sqrt(2)*hbar*k0**2)*(1-ep)**1.5
    
    return np.exp(-0.5*(termino1*termino2+termino3))

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


#definiendo energias adimensionales
Ne=100
E=np.linspace(0,2,Ne)

eta=[1E-7,1E-5,1E-7]# eV/Fermi

#Constantes usadas


Vo=1.8E-6
m=0.5E6# MeV
hbar=197.327# MeV*fermi hbar*c



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
a=20.8E5#fermi
Nx=1000

                

for i in range(len(E)):
    
    #defining the array for the intergration
    
    tpoint=np.sqrt(2*m*epsilon[i]*Vo*(E[i]-1))/eta[0]
    x=np.linspace(tpoint,a,Nx)
    dx=abs(x[0]-x[1])
    
    #suma
    suma=0
    
    expo=function(m,eta[0],Vo,E[i],hbar,epsilon[i])
        
    Kfunc=k(m,eta[0],Vo,E[i],hbar,epsilon[i])
    
    #se empieza la integral
    for j in range(len(x)):
        
        suma+=dx*(1/Kfunc[j])*expo[j]
    
    #este es el tiempo hallado
    T[i]=suma


        

        













