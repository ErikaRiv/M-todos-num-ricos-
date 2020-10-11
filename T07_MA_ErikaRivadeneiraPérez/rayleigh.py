#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:38:34 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%%----------------------RAYLEIGH---------------------------------------------
#A : matriz a la que se desea encontrar el autovalor más grande
#max_iter : Número máximo de iteraciones aceptadas 
#tol : Tolerancia aceptada

#AUTOVALOR MAS CERCANO A UN NUMERO INICIAL 
def rayleigh(A, max_iter, tol): 
    n,m = np.shape(A)
    I = np.eye(n)
    vo = np.zeros(n)
    vo[n-1] = 1
    vo = vo / np.linalg.norm(vo)
    cont = 0
    lambd_aprox = 1
    while( cont < max_iter):
        v1 = np.linalg.solve(A - lambd_aprox*I,vo)
        v1 = v1/ np.linalg.norm(v1)
        lambd = v1@A@v1
        error = abs(lambd_aprox-lambd)/abs(lambd)
        lambd_aprox=lambd
        vo[:] = v1
        cont +=1
        if (error < tol):
            print("Autovalor: ", lambd_aprox)
            print("\nAutovector: ", vo)
            print("\nIteraciones: ", cont)
            print("\nError: ", error)
            return lambd_aprox, vo, cont, error
    
#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A1=np.loadtxt("Eigen_50x50.txt", skiprows=1)#A50x50
#%% TOLERANCIA E ITERACIONES SEGURIDAS
tol = 1e-7
max_iter=10000
#%%------------------------RAYLEIGH EN MATRIZ 3x3------------------------------
print("------------RAYLEIGH EN MATRIZ 3x3-------------\n")
rayleigh(A,max_iter, tol)

#%%----------------------RAYLEIGH  EN MATRIZ 50x50-----------------------------
print("------------RAYLEIGH EN MATRIZ 50x50-------------\n")
rayleigh(A1,max_iter, tol)

    