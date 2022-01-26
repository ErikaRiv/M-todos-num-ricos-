#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 23:34:46 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%% METODO DE FACTORIZACION QR
def QR(A):
    m, n = np.shape(A)
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    A = np.transpose(A)
    R[0][0]=np.linalg.norm(A[0], 2)
    Q[:,0]=A[0]/R[0][0]
    lambds = []
    for i in range(1,n):
        Q[:,i]=A[i]
        for j in range(0,i):
            R[j][i]=np.dot(Q[:,j], Q[:,i])
            Q[:,i]=Q[:,i]-(R[j][i] * Q[:,j])
        R[i][i]=np.linalg.norm(Q[:,i], 2)
        Q[:,i]=Q[:,i]/R[i][i]
    for i in range(n):
        lambds.append(R[i,i])
    
    return Q, R
#%% SUSTITUCION HACIA ATRAS
def solve_triang_superior( matriz, b ):
    n=len(matriz) #Número de elementos de A
    x = np.zeros(n)#Creamos vector de ceros para guardar soluciones
    for i in range(n-1,-1,-1):#recorremos la matriz desde el final hasta el inicio
        sumj = 0 #inicializamos las sumas
        for j in range(i+1,n): 
            sumj += matriz[i,j]*x[j] #Realizamos las sumas de cada iteracion
        x[i] = (b[i]-sumj)/matriz[i,i]
    return x #Retornamos las raices 

#%% SOLVER QR
#Rx=QTb.
def solver_QR(A,b):   
    Q=QR(A)[0]
    R=QR(A)[1]
    b_new = np.transpose(Q)@b
    x=solve_triang_superior( R, b_new )
    error = np.linalg.norm(A@x-b)
    print("Solución: \n", x)
    print("Error: ", error)
    return x, error
#%%QR EIGEN VALUES
#A matriz para encontrar autovalores
#tol: Tolerancia o presición aceptada
#max_iter: número máximo de iteraciones
def values_QR(A, tol, max_iter):
    n,m = np.shape(A)
    Vect = np.eye(n)
    cont = 0
    while(cont < max_iter):
        Q,R = QR(A)
        A1 = R@Q
        Vect = np.transpose(Q)@Vect
        error = np.linalg.norm(A-A1)
        if (error<tol):
            return np.diag(A1),Vect,cont,error
        A = A1.copy()
        cont+=1
    return np.diag(A1),Vect,cont,error
#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A1=np.loadtxt("Eigen_50x50.txt", skiprows=1)#A50x50
A2=np.loadtxt("M_sys_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A3=np.loadtxt("M_sys_125x125.txt", skiprows=1)#A50x50
b2=np.loadtxt("V_sys_3x1.txt", skiprows=1)#A50x50
b3=np.loadtxt("V_sys_125x1.txt", skiprows=1)#A50x50
#%%
tol = 1e-10
max_iter = 10000
#%%#IMPRIMIENDO RESULTADOS DE FACTORIZACION
print("------Factorización QR de la matriz A3X3------\n")
print("Q = ",QR(A)[0])
print("\nR = ",QR(A)[1])
vasl2= values_QR(A,tol, max_iter)
print("\nAutovalores = ", vasl2[0], "\nIteraciones: ",vasl2[2],"\nError: ", vasl2[3])

#%%
print("------Factorización QR de la matriz A50X50------\n")
print("Q = ",QR(A1)[0])
print("\nR = ",QR(A1)[1])
vals =  values_QR(A1,tol, max_iter)
print("\nAutovalores = ", vals[0], "\nIteraciones: ",vals[2],"\nError: ", vals[3])
#%%IMPRIMIENOO SOLUCION DE SISTEMAS Ax=b
print("------Solución del sistema A3x3X=b3x1------\n")
solver_QR(A,b2)
#%%
print("------Solución del sistema A125x125X=b125x1------\n")
solver_QR(A3,b3)
