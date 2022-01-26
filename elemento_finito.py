#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 23:47:50 2020

@author: erikarivadeneira
"""

#LIBRERIAS
import numpy as np
import matplotlib.pylab as plt
#%% ELEMENTO FINITO 
#x: vector de puntos x
#y: vector de imagenes de x
#xo: parámetro de penalización
def elemento_finito(x,y,xo):
    polinomio = np.zeros(len(xo))
    n=len(x)
    cont=0
    for i in range(n-2):
        while(xo[cont]<=x[i+1]): 
            n1 = (xo[cont]-x[i+1])/(x[i]-x[i+1])
            n2 = (xo[cont]-x[i])/(x[i+1]-x[i])
            yj = y[i]*n1+y[i+1]*n2
            polinomio[cont]=yj
            cont+=1
    m=len(polinomio)
    polinomio[m-1]=y[n-1]
    plt.plot(x,y, '*', color= 'indigo')        
    plt.plot(xo,polinomio, ls='-',color = 'lightslategray')
    plt.show()

#%% CARGANDO DATOS
data=np.loadtxt("datos.txt",delimiter= ',')
x=data[:,0]
y=data[:,1]
x_min = min(x)
x_max = max(x)
#%%PARAMETRO DE PENALIZACION = 10
xo = np.linspace(x_min,x_max,10)
elemento_finito(x,y,xo)
#%%PARAMETRO DE PENALIZACION = 40
xo1 = np.linspace(x_min,x_max,40)
elemento_finito(x,y,xo1)
#%%PARAMETRO DE PENALIZACION = 80
xo2 = np.linspace(x_min,x_max,80)
elemento_finito(x,y,xo2)