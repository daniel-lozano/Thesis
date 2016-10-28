import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def func(x,a):
    return np.sin(x)+a


x=np.linspace(0,np.pi,1000)
y=np.sin(x)+1

suma=0

for i in range(len(x)):
    suma+=(x[1]-x[0])*y[i]
print("resultado usando suma tipica",suma)

funcion= lambda x: func(x,1)

resultados=integrate.quad(funcion,0,np.pi,limit=100)

print("resultado usando cuadratura gaussiana",resultados[0], " Error de Quad ", resultados[1])