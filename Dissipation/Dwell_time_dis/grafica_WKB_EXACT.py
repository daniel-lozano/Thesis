import numpy as np
import matplotlib.pyplot as plt

E=np.linspace(0.1,0.99)
T_wkb=np.zeros(len(E))
T_exact=np.zeros(len(E))
#definiendo las constantes

a=20.8E5#fermi
Vo=1.8E-6
m=0.5E6# MeV
hbar=197.327# MeV*fermi =[hbar*c]
eta=0

#Llenando T_wkb

for i in range(len(E)):
    k=np.sqrt(2*m*Vo*(1-E[i]))/hbar
    result=(m/(2*hbar*pow(k,2)))*( 1-np.exp(-2*k*a-eta*pow(a,2)/hbar)*(2/(1+eta*a/(hbar*k))-1) )
    T_wkb[i]=result*6.58*pow(10,-16)/hbar

#Llenando T_wkb

for i in range(len(E)):
    k=np.sqrt(2*m*Vo*(1-E[i]))/hbar
    kp=np.sqrt(2*m*Vo*E[i])/hbar

    result=2*m*(a*k*( pow(kp,2)-pow(k,2) ) + np.sinh(2*kp*a)*((pow(kp,2)+pow(k,2))/2*k*kp)   )/(( np.cosh(kp*a)*(pow(kp,2)+pow(k,2) ))**2 - (pow(kp,2)-pow(k,2) )**2  )
    T_exact[i]=result






plt.plot(E,T_wkb,label="$ WKB $")
plt.xlabel("$ E/V_0 $")
plt.ylabel("$ Dwell\ Time\ $")
plt.legend(loc=2)
plt.show()

