#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:43:51 2020

@author: erikarivadeneira
"""
#Librerias

import numpy as np

#%%
#-----------------------------GAUSS-SEIDEL--------------------------------------
#A: Matriz a resolver
#b: Vector independiente
#tol: Tolerancia que el usuario desee
#max_ire: Número máximo de iteraciones
#El vector inicial x que se considera es un vector de ceros
def Gauss_seidel(A,b,tol,max_iter):
    error = 1
    m,n = np.shape(A)
    x = np.zeros((m))
    cont = 0
    #comp = []
    #comp=np.zeros((n))
    for i in range(m):
        if(A[i,i]!=0):
            while(cont<max_iter and error>tol):
                suma=0
                cont+=1
                for i in range(m):
                    suma = 0
                    for j in range(n):
                        if (j != i):
                            suma+=A[i,j]*x[j] 
                    x[i]=(b[i]-suma)/A[i,i]
                   # print("x[",i,"]:",x[i])
                print("Iteraciones:  ", cont)
                error = np.linalg.norm(np.dot(A,x)-np.transpose(b))
    print("Error: ", error)
    print("\nEl programa terminó en la iteranción: ", cont)
    print("\nSolución: ")
    print(x)
    return x
    
#%%
#--------------------------------JACOBI----------------------------------------
#A: Matriz a resolver
#b: Vector independiente
#tol: Tolerancia que el usuario desee
#max_ire: Número máximo de iteraciones
#El vector inicial x que se considera es un vector de ceros
def Jacobi(A,b,tol,max_iter):
    m,n = np.shape(A)
    x = np.ones((n))
    xo = np.zeros((n))
    cont = 1
    error=1
    #comp=np.zeros((n))
    for i in range(m):
        if(A[i,i]!=0):
            while(cont<max_iter and error>tol):
                cont=cont+1
                for k in range(m): 
                    suma=0
                    for j in range(m):
                        if (k != j):
                            suma+=A[k,j]*xo[j] 
                    x[k]=(b[k]-suma)/A[k,k]
                    #print("x[",k,"]:",x[k])
                print("Iteraciones:  ", cont)
                xo[:] = x
                error = np.linalg.norm(np.dot(A,xo)-np.transpose(b))
    print("El programa terminó en la iteranción: ", cont)
    print("\nError: ", error)
    print("\nSolución: ")
    print(x)
    return x
    
#%%
#------------------------------RESULTADOS--------------------------------------
#Probando métodos iterativos para solución de sistemas lineales
#Cargando matrices de dimensión n=3
A=np.loadtxt("M_sys_3x3.txt",delimiter= ' ', skiprows=1)
b=np.loadtxt("V_sys_3x1.txt",delimiter= ' ', skiprows=1)
tol=0.00001
max_iter=1000
#Cargando matrices de dimensión n=125
A1=np.loadtxt("M_sys_125x125.txt",delimiter= ' ', skiprows=1)
b1=np.loadtxt("V_sys_125x1.txt",delimiter= ' ', skiprows=1)

#Sugerencia de tolerancia y máximo número de iteraciones 
tol=0.0000001
max_iter=1000
#%%--------------PRUEBA CON GAUSS-SEIDEL PARA A3X3 Y b3X1---------------------
print("Gauss-Seidel\n")
Gauss_seidel(A,b,tol,max_iter)

#%%-----------------PRUEBA CON JACOBI PARA A3X3 Y b3X1------------------------
print("\nJacobi\n")

Jacobi(A,b,tol,max_iter)

#%%----------PRUEBA CON GAUSS-SEIDEL PARA A125X125 Y b125X1--------------------
print("Gauss-Seidel\n")
Gauss_seidel(A1,b1,tol,max_iter)
#%%-------------PRUEBA CON JACOBI PARA A125X125 Y b125X1----------------------
print("\nJacobi\n")
Jacobi(A1,b1,tol,max_iter)

