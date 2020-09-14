#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 21:36:54 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from scipy import sparse

#%%
#----------------------------------CHOLESKY------------------------------------  
def cholesky(A):
    #Primero comprobamos que la matriz sea simétric y def positiva
    if (es_simetrica(A)==True & es_def_pos(A)==True):
        n = len(A)
        for i in range(n):
            suma = A[i,i]
            for k in range(i):
                suma -= A[k,i]**2
            A[i,i] = np.sqrt(suma)
            for j in range(i+1, n):
                suma = A[i,j]
                for k in range(i):
                    suma -= A[k,i]*A[k,j]
                A[i,j]= suma / A[i,i]
    
       	for j in range(n):
            for i in range(n):
                if(i > j):
                    A[i,j] = 0.0
        Lt = np.transpose(A)
        L = A
        print("La factorización LLT de A es:\n L=",Lt, "\n\nLT=",L)
        return Lt, L, Lt*L
        
        
    else:
        print("La matriz no es definida positiva o simétrica")
        
#%% 
#-----------------------------FACTORIZACION LDLT-------------------------------  
def ldlt(A):
    #Primero comprobamos que la matriz sea simétrica
    if( es_simetrica(A) == True):
        n,m = np.shape(A)
        L = np.eye(n)
        D = np.zeros((n, 1))
        for i in range(n):
            D[i] = A[i, i] - np.dot(L[i, 0:i] ** 2, D[0:i])
            for j in range(i + 1, n):
                L[j, i] = (A[j, i] - np.dot(L[j, 0:i] * L[i, 0:i], D[0:i])) / D[i]
        D = np.eye(n) * D
        print("La factorización LDLT de A es:\n L=",L, "\nD=",D,"\nLT=",np.transpose(L))
        return L, D, np.transpose(L)
    else:
        print("La matriz no es simétrica")
#%%
#-------------------------FUNCIONES ADICIONALES--------------------------------
#Verifica que la matriz sea simétrica
def es_simetrica(A):
    sim = "T"
    n,m = np.shape(A)
    for i in range(n):
        for j in range(m):
            if A[i,j]!=A[j,i]:
                sim = "F"
    if(sim == "T"):
         return True
    else:
        return False
#Función que verifica que la matriz sea definida positiva
def es_def_pos(A):
       return np.all(np.linalg.eigvals(A) > 0)
   
#Creando matriz A de la ecuación del calor 
def CrearMatrizCalor(n):
    dl = -1 * np.ones(n-1)
    du = -1 * np.ones(n-1)
    d0 = 2 * np.ones(n-1)
    d = np.vstack((dl, d0, du))
    A = sparse.spdiags(d, (-1, 0, 1), n-1 , n-1 )
    A = A.todense()
    return A

    
 #%%   
#----------PROBANDO CHOLESKY Y LDLT PARA N=4-----------------------------------
#Los resultados de matrices grandes se encuentran en archivos txt
a=ldlt(CrearMatrizCalor(4))
#Para comprobar que el resultado esté bien hecho hacemos 
print("Verificación:\nA=LDLT=",a[0]@a[1]@a[2])
print("\n")
print("Factorización de Cholesky de A es:\n")
b=cholesky(CrearMatrizCalor(4))
print("Verificación:\nA=LLT=",b[0]@b[1])

#Escribir archivo txt para almacenar resultados de matrices grandes
#c=cholesky(CrearMatrizCalor(100))
#np.savetxt("L100_CHOLESKY.txt",c[0],fmt="%f"	)
    
