import numpy as np
import matplotlib.pyplot as plt



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def k(m,eta,Vo,E,hbar,x):#sheck
    
    k0=np.sqrt(2*m*Vo*(1-E))/hbar
    
    return np.sqrt(k0**2+(2*eta*k0*x)/hbar +(eta*x/hbar)**2 )


def function(m,eta,Vo,E,hbar,x):
    
    k0=np.sqrt(2*m*Vo*(1-E))/hbar
    
    termino1=-0.5*(2*x+2*hbar*k0/eta)
    
    termino2=np.sqrt(k0**2 +(2*eta*k0*x)/hbar +(eta*x/hbar)**2  )
    
    termino3= (1/eta)*(4*hbar*k0**2)
    
    return np.exp(termino1*termino2+termino3)

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


#definiendo energias adimensionales
Ne=100
E=np.linspace(0,0.99,Ne)

#Constantes usadas

a=20.8E5#fermi
Vo=1.8E-6
m=0.5E6# MeV
hbar=197.327# MeV*fermi =[hbar*c]
print("\nmass=",m, " MeV")
print("L=",a," fermi")
print("Vo=",Vo, " MeV")
print("hbar=",hbar," Mev*fermi\n")


'''
Tenemos que tener que eta*l << hbar*k, para satisfacerlo tomaremos el k mas pequeño posible correspondiente a E=0
'''

hbar_k_l=np.sqrt(2*m*Vo*(1-E[-1]))/a

print("Condition:\neta <<",hbar_k_l)

eta=[(hbar_k_l)*10**-2,hbar_k_l*10**-3,hbar_k_l*10**-4] # eV/Fermi

print("eta=",eta)


T=np.zeros(len(E))

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comenzando el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


Nx=1000

print("\nEnergies used\n",E)
k=np.zeros(len(E))

for i in range(len(E)):
    k[i]=np.sqrt(2*m*Vo*(1-E[i]))/hbar

plt.plot(E,k)
plt.show()


#defining the array for the intergration

x=np.linspace(0,a,Nx)
dx=abs(x[0]-x[1])


print("probando",function(m,eta[0],Vo,E[1],hbar,x[0]))


for i in range(len(E)):
    
    #suma
    suma=0
    
    #se empieza la integral
    
    
    for j in range(len(x)):
        
        expo=function(m,eta[0],Vo,E[i],hbar,x[j])
    
        Kfunc=k(m,eta[0],Vo,E[i],hbar,x[j])

        suma+=dx*(1/Kfunc)*expo
    
    #este es el tiempo hallado
    T[i]=suma*m/hbar

        

        













