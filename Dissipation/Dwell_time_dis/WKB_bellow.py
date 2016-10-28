import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiempo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def integrand(x,k,eta,hbar,m,a):
    exp=np.exp(-2*k*x*(1+eta*x/2))
    div=1+eta*x
    return exp/div

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


#definiendo energias adimensionales
Ne=100
E=np.linspace(0.01,0.99,Ne)

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

eta=np.linspace(0,0.5,3)*10**(-7)/hbar_k_l # eV/Fermi

print("eta=",eta)


T=np.zeros(len(E))

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comenzando el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

#defining the array for the intergration


for j in range(len(eta)):

    for i in range(len(E)):
    
        k=np.sqrt(2*m*Vo*(1-E[i]))/hbar
        kp=np.sqrt(2*m*Vo*(E[i]))/hbar
        
        N_factor=1#np.sqrt(k*kp/(k**2+kp**2))
    
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
    
        resultados=integrate.quad(funcion,0,1,limit=100)

        T[i]=(N_factor**2)*(6.58*10**(-16))*(m*a/(hbar*k))*(resultados[0])/hbar

        
    plt.plot(E,T)

#plt.figure(figsize=(20,10))
plt.show()
plt.savefig("Integral_disp_vel.png")










