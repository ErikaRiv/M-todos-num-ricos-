#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 21:25:06 2020

@author: erikarivadeneira
"""

import numpy as np
import matplotlib.pyplot as plt
#%% FUNCION PARA SPLINES CUBICOS
#x: vector de puntos del eje x
#y: vector de inmágenes de x
#puntos: # de puntos de entrenamient0
#funcion: funcion real
def splines_cubicos(x,y,puntos,funcion):
    num_puntos = puntos
    n = len(x)
    Hx = np.diff(x)
    Hy = np.diff(y)
    Hn = n - 1
    Alfa = 1 / Hx[1 : Hn - 1]
    Gamma = 1 / Hx[1 : Hn - 1]
    Beta = 2 * (1 / Hx[:Hn - 1] + 1 / Hx[1:])
    dF = Hy / Hx
    Delta = 3 * (dF[1:] / Hx[1:] + dF[:Hn-1] / Hx[:Hn-1])
    TDM = np.diag(Alfa, k=-1) + np.diag(Beta, 0) + np.diag(Gamma, +1)
    B = np.linalg.solve(TDM, Delta)
    B = np.hstack([0, B, 0])    
    C = (3*dF - B[1:] - 2 * B[:Hn]) / Hx
    D = (B[:Hn] + B[1:] - 2 * dF) / (Hx ** 2)
    x_paso = (x[n-1] - x[0]) / num_puntos
    x_puntos = []
    x_base = x[0]
    for i in range(num_puntos):
        x_puntos.append(x_base+x_paso*i)
    y_puntos = []
    for val in x_puntos:
        for i in range(n-1):
            if ((val >= x[i]) and (val <= x[i+1])):
                y_val = y[i] + B[i] * (val - x[i]) + C[i] * ((val - x[i]) ** 2) + D[i] * ((val - x[i]) ** 3)
                y_puntos.append(y_val)
    
    plt.plot(np.arange(x[0],x[n-1],0.00001),funcion,color="chartreuse")
    plt.plot(x_puntos, y_puntos,"-g")
    plt.plot( x_puntos, y_puntos, "*",color="rebeccapurple")
    plt.show()
   
#%%PROBANDO INTERPOLACION POR SPLINES EN FUNCIONES DADAS
print("Función f(x)=sin(x)")
x = np.linspace(-np.pi,np.pi,15)
y = np.sin(x) 
print("\n -----15 puntos de prueba-------")
splines_cubicos(x,y,15,np.sin(np.arange(x[0],x[len(x)-1],0.00001)))
#%%
print("\n -----150 puntos de prueba-------")
splines_cubicos(x,y,150,np.sin(np.arange(x[0],x[len(x)-1],0.00001)))
#%%FUNCION 2 
print("\nFunción f(x)=x^2(5x-3)-2x^4+4x-5")
x2 = np.linspace(-2,4,14)
y2 = x2**2*(5*x2-3)-2*x2**4+4*x2-5
print("\n -----15 puntos de prueba-------")
xo = np.arange(x2[0],x2[len(x2)-1],0.00001)
funcion = xo**2*(5*xo-3)-2*xo**4+4*xo-5
splines_cubicos(x2,y2,15,funcion)
#%%
print("\n -----150 puntos de prueba-------")
splines_cubicos(x2,y2,150,funcion)
#%%FUNCION 3 
print("\nFunción f(x)=e^-x^2")
x3 = np.linspace(-1,1,14)
y3 = np.exp(-x3**2)
print("\n -----15 puntos de prueba-------")
xo3 = np.arange(x3[0],x3[len(x3)-1],0.00001)
funcion2 = np.exp(-xo3**2)
splines_cubicos(x3,y3,15,funcion2)
#%%
print("\n -----150 puntos de prueba-------")
splines_cubicos(x3,y3,150,funcion2)