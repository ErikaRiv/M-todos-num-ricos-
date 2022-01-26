#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 21:06:49 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
import sympy as sym

#%% NEWTON ANIDADO 
#a: vector de coeficientes 
#xi:vector de centros
#c: centro a evaluar en el polinomio de Newton anidado 
def evaluacion_polinomio_Newton(a,xi,c):
    v=0
    n=len(a)
    v = a[n-1]
    for i in range(n-2,-1,-1):
        v = v*(c-xi[i])+a[i]
    print("Evaluación del polinomio de Newton en el centro ",c)
    print("P",n-1,"(",c,") = ", v)
    return v   

#%%ESCRIBIENDO EN LA FORMA DEL POLINOMIO 
#a: vector de coeficientes 
#xi:vector de centros
def forma_pol_Newton(a,xi):
    x = sym.Symbol('x')
    n=len(a)
    v = a[n-1]
    polinomio = 0
    for i in range(n-2,-1,-1):
        v = v*(x-xi[i])+a[i]
    polinomio += v
    polinomio_simple = sym.expand(polinomio)
    print("\nPolinomio en forma de Newton anidada: \n", polinomio)
    print("\nPolinomio en forma simple: \n",polinomio_simple)

#%% Probando funciones con el ejemplo 8.1 
a=np.array([0,1,2,3,4])
xi = np.array([0,1,2,3])

evaluacion_polinomio_Newton(a,xi,4)
forma_pol_Newton(a,xi)

np.array([[1, 2, 3, 4],
          [7, 7, 7, 7],
          [7, 7, 7, 7]])
