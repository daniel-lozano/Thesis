import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiempo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def integrand(x,k,eta,hbar,m,a):
    
    factor1= - 2*k*a*x

    factor2= - (eta/hbar)*pow(x*a,2)

    exp= np.exp(factor1+factor2)

    div= hbar*k + eta*x*a

    return a*m*exp/div

'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


#definiendo energias adimensionales
Ne=100
E=np.linspace(0.01,0.99,Ne)

#Constantes usadas

a=20.7989E5  #fermi
Vo=1.8E-6 # MeV
m=0.51     # MeV
hbar=197.327# MeV*fermi =[hbar*c]

print("\nmass=",m, " MeV")
print("L=",a," fermi")
print("Vo=",Vo, " MeV")
print("hbar=",hbar," Mev*fermi\n")


'''
Tenemos que tener que eta*l << hbar*k, para satisfacerlo tomaremos el k mas pequeño posible correspondiente a E=0
'''


eta=np.linspace(0,5E-11,2) # eV/Fermi

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
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
    
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
    
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)

        T[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)
        
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV $"

        
    plt.plot(E,T,label=nombre)

for i in range(len(T)):
    
    print(E[i],T[i])






#plt.figure(figsize=(20,10))

plt.legend(loc=2)
plt.xlabel("$  E/V_0 $",size=20)
plt.ylim(0,1E-15)
plt.ylabel(" $ Time $ ",size=20)
plt.savefig("Integral_disp_vel.png")

plt.show()









