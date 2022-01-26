#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 19:23:33 2020

@author: erikarivadeneira
"""

#LIBRERIAS
import numpy as np
import matplotlib.pylab as plt
#%%AJUSTE POR MINIMOS CUADRADOS
#x: vector de puntos x
#y: vector de imagenes de x
#tipo: lineal ->1, parabola -> 2, exponencial ->3
def min_squares(x,y,tipo):
    n = len(x)
    x_min = min(x)
    x_max = max(x)
    xo = np.linspace(x_min,x_max,n)
    if tipo == 1: #AJUSTE LINEAL
        A = np.ones((len(x),2))
        A[:,0]=x
        B=y
        ATb = np.dot(np.transpose(A),B)
        ATA = np.linalg.inv(np.dot(np.transpose(A),A))
        X=np.dot(ATA,ATb)
        a=X[0]
        b=X[1]
        yo = a*xo + b
    if tipo == 2: #AJUSTE DE PARABOLA
        A = np.ones((n,3))
        A[:,0] = x**2
        A[:,1] = x
        B = y
        ATb = np.dot(np.transpose(A),B)
        ATA = np.linalg.inv(np.dot(np.transpose(A),A))
        X=np.dot(ATA,ATb)
        a=X[0]
        b=X[1]
        c=X[2]
        yo=a*xo**2+b*xo+c
        
    if tipo == 3: #AJUSTE EXPONENCIAL
        A = np.ones((len(x),2))
        A[:,0]=x
        B=y
        ATb = np.dot(np.transpose(A),B)
        ATA = np.linalg.inv(np.dot(np.transpose(A),A))
        X=np.dot(ATA,ATb)
        a=X[0]
        b=X[1]
        yo=a+b*np.exp(x)
    plt.plot(x,y,'*',color= 'indigo')
    plt.plot(xo,yo,color='lightslategray')
    plt.show()
    yo =yo/np.linalg.norm(yo)
    y_real = 2.15*np.exp(0.758*x)
    y_real = y_real/np.linalg.norm(y_real)
    error = (np.square(yo-y_real)).mean()#Cálculo de error
    
    return error
        
#%% CARGANDO DATOS
data=np.loadtxt("datos.txt",delimiter= ',')
y=data[:,1]
x=data[:,0]

#%%AJUSTE LINEAL
error=min_squares(x,y,1)
print("Error de ajuste lineal: ", error)
#%%#AJUSTE PARABÓLICO
error2=min_squares(x,y,2)
print("Error de ajuste parabólico: ", error2)

#%%AJUSTE EXPONENCIAL
error3=min_squares(x,y,3)
print("Error de ajuste exponencial: ", error3)

#%%

