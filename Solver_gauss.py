#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:59:01 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from Matriz_diagonal import crearMatriz, imprimir_sol
from Solver_triangular_sup import solve_triang_superior
#%%

#Creamos función para eliminacion gaussiana sin pivoteo 
def eliminacion_gaussiana(A, b):
    matrix = np.insert(A,A.shape[1],np.transpose(b),1) #Concateno A y b para operar simultaneamente
    np.asarray(matrix) #aseguro el arreglo
    matrix = matrix.astype(float) #hago flotantes a los datos del arreglo
    n,m = matrix.shape
    #empezamos la eliminación
    for i in range(0,n):#row
        for j in range(i+1,n):
            if matrix[j,i] != 0.0:
                #operamos
                multiplier = matrix[j,i]/matrix[i,i]
                matrix[j,i:m]=matrix[j,i:m] - multiplier*matrix[i,i:m] 
    #Resolvemos Ax=b
    A = matrix[:,:n]
    b = matrix[:,n:n+1]
    return solve_triang_superior(A,b)

#%%
#--------------Para eliminación Gaussiana---------------------
#----------sin pivoteo-------------------------------------------
print("\n---------------Eliminación Gaussiana sin pivoteo-------------------\n")
#Para matriz pequeña
A_small = crearMatriz("M_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_small = np.matrix(A_small)#transformamos el array en numpy matriz

b_small = crearMatriz("V_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_small = np.matrix(b_small)#transformamos el array en numpy matriz

Sol_small = eliminacion_gaussiana(A_small,b_small) #Encontramos la solucion del sistema Ax=b para la matriz
imprimir_sol( Sol_small, A_small,b_small)#Imprimimos solucion
#%%
#Para matriz grande
A_large = crearMatriz("M_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_large = np.matrix(A_large)#transformamos el array en numpy matriz

b_large = crearMatriz("V_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_large = np.matrix(b_large)#transformamos el array en numpy matriz

Sol_large = eliminacion_gaussiana(A_large,b_large) #Encontramos la solucion del sistema Ax=b para la matriz
imprimir_sol( Sol_large, A_large,b_large)#Imprimimos solucion