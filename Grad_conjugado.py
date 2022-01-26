#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 23:33:15 2020

@author: erikarivadeneira
"""
#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%%--SOLUCION DE SISTEMA DE ECUACIONES POR EL METODO DE GRADIENTE CONJUGADO----
#DATOS DE ENTRADA PARA EL SISTEMA Ax=b:
#A = Matriz A del sistema
#b = vector de términos independientes del sistema
#tol = Tolerancia o presición
#max_iter = Número máximo de iteraciones
def Grad_conjugado(A, b, tol, max_iter):
    n = len(b)
    x0 = np.zeros(n)
    r= np.dot(A, x0) - b        
    d =  -r                     
    cont = 0
    error = 1
    while(cont < max_iter and error > tol ):
        alpha = np.dot(np.transpose(r), r) / np.dot(np.dot(np.transpose(d), A), d)
        x0 = x0 + d*(alpha)
        r1 = r + np.dot(A*(alpha), d)
        error = np.linalg.norm(r1)
        betai = np.dot(np.transpose(r1), r1) / np.dot(np.transpose(r), r)
        d = - r1 + (betai)*d
        r = r1
        cont+=1
    print("Solución:\n ", x0)
    print("\nIteraciones: ", cont)
    print("\nError:", error )
    return x0, cont, error

#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
b=np.loadtxt("V_sys_3x1.txt", skiprows=1)#b3x1
A1=np.loadtxt("M_sys_125x125.txt", skiprows=1)#A125x125
b1=np.loadtxt("V_sys_125x1.txt", skiprows=1)#b125x1

#MAXIMO NUMERO DE ITERACIONES Y TOLERANCIA SUGERIDA 
max_iter=10000
tol= 1e-10
#%%PROBANDO METODO DE GRADIENTE CONJUGADO 
print("--------GRADIENTE CONJUGADO-------------")
print("\n----Solucion Ax=b para matriz A de dimension 3x3------\n")
Grad_conjugado(A, b, tol, max_iter)
#%%
print("\n----Solucion Ax=b para matriz A de dimension 125x125------\n")
Grad_conjugado(A1, b1, tol, max_iter)