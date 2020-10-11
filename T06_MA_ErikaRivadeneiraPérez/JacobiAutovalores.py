#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 00:30:58 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np

#%%----------------------------METODO DE JACOBI--------------------------------
#A : matriz a la que se desea encontrar el autovalor más grande
#max_iter : Número máximo de iteraciones aceptadas 
#tol : Tolerancia aceptada
def Jacobi_autovalores(A,max_iter,tol):
    n,m = np.shape(A)
    cont = 0
    error = 1
    vectores = np.eye(n)
    while(cont<max_iter):
        R = np.eye(n)
        L = np.tril(A)
        U = A-L
        Arg_max = np.amax(abs(U))
        (i,j) = np.where(abs(U)==Arg_max)
        i = i[0]
        j = j[0]
        error = Arg_max
        suma = 0 
        for k in range(n):
            suma += A[k,k]*A[k,k]
        if (error < tol):
            break
        if(A[i,i]==A[j,j]): #Encontramos el ángulo theta
            theta = np.pi/4
        else:
            theta = 0.5*np.arctan(2.0*A[i,j]/((A[i,i]-A[j,j])))             
        R[i,i] = np.cos(theta) #Empezamos rotaciones
        R[j,j] = np.cos(theta)    
        R[i,j] = np.sin(theta)
        R[j,i] = -np.sin(theta)
        vectores = np.dot(vectores,np.transpose(R)) 
        prod = np.dot(R,A)
        A = np.dot(prod, np.transpose(R))
        cont+=1
    valores = np.diag(A)
    print("Autovalores: ", valores)
    print("\n Autovectores: ", vectores)
    print("Iteración: ", cont)
    print("Error: ", error)
    return valores, vectores
                
#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A1=np.loadtxt("Eigen_50x50.txt", skiprows=1)#A50x50
A2=np.loadtxt("Eigen_500x500.txt", skiprows=1)#A500x500
A3=np.loadtxt("Eigen_1000x1000.txt", skiprows=1)#A1000X1000

#%% TOLERANCIA E ITERACIONES SEGURIDAS
tol = 1e-7
max_iter=10000
#%%--------------------------JACOBI EN MATRIZ 3x3------------------------------
print("------------JACOBI EN MATRIZ 3x3-------------\n")
Jacobi_autovalores(A,max_iter, tol)

#%%------------------------JACOBI EN MATRIZ 50x50------------------------------
print("------------JACOBI EN MATRIZ 50x50-------------\n")
Jacobi_autovalores(A1,max_iter, tol)
