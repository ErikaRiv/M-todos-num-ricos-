#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 18:40:18 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%%----------------------------METODO DE JACOBI--------------------------------
#A : matriz a la que se desea encontrar autovalores
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
    return valores, vectores
#%%----------------------ITERACION EN SUBESPACIOS-----------------------------
#A : matriz a la que se desea encontrar autovalores
#R : número de autovalores deseados 
#max_iter : Número máximo de iteraciones aceptadas 
#tol : Tolerancia aceptada
def iteracion_subespacio(A, R, max_iter, tol ):
    n,m = np.shape(A)
    cont = 0
    Phi = np.eye(n,R) # Inicializo matriz de subespacio
    while( cont < max_iter):
        for i in range(R):
            vo = Phi[:,i].copy()
            aux = vo.copy()
            for j in range(i):
                aux -= (Phi[:,j]@vo)*Phi[:,j]
            vo[:]=aux
            vo = vo / np.linalg.norm(vo)
            v1 = np.linalg.solve(A,vo)
            v1 = v1 / np.linalg.norm(v1)
            v1_transp = np.transpose( v1 )
            Phi[:,i]= v1_transp
        Q = np.transpose(Phi)@A@Phi
        error = 0
        cont +=1
        for i in range(R):
            error +=sum(abs(Q[i,(i+1):R]))
        if (error < tol):
            lambds = np.diag(Q)
            vects = Phi
            break
        val, vect = Jacobi_autovalores(Q,max_iter,tol) #Retorno los autovlores y autovectores de interés
        Phi = Phi@np.transpose(vect)
        
    print("---------ITERACIÓN EN SUBESPACIO---------\n")
    print("Autovalores:  \n", lambds)
    print("\nAutovectores:  \n", vects)
    print("\nError:  ", error)
    print("\nIteraciones: ", cont)

    return lambds, vects 
        
#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A1=np.loadtxt("Eigen_50x50.txt", skiprows=1)#A50x50
#%% TOLERANCIA E ITERACIONES SEGURIDAS
tol = 1e-7
max_iter=10000
#%%----------------ITERACION EN SUBESPACIO EN MATRIZ 3x3-----------------------
print("-----------ITERACION EN SUBESPACIO EN MATRIZ 3x3-------------\n")
iteracion_subespacio(A, 2, max_iter, tol )
#%%----------------ITERACION EN SUBESPACIO EN MATRIZ 50x-----------------------
print("------------ITERACION EN SUBESPACIO 50x50-------------\n")
iteracion_subespacio(A1, 10, max_iter, tol )
#np.savetxt("Autovals_50x50.txt",R1[0], fmt="%f"	)





