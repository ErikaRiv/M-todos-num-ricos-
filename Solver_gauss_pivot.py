#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 20:01:14 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from Matriz_diagonal import crearMatriz
from Solver_triangular_sup import solve_triang_superior
#%%
#Creamos función para eliminación gaussiana con pivoteo
def eliminacion_gaussiana_piv(A, b):
    matrix = np.insert(A,A.shape[1],np.transpose(b),1)
    n,m = np.shape(matrix) #numero de filas

	# Create an upper triangular matrix
    for i in range(0,n): # for i in [0,1,2,..,n]
		#Encontramos el valor más grande en las columnas 
        maxElem = abs(matrix[i,i])
        maxRow = i
        for k in range(i+1, n):
            if( abs(matrix[k,i])>maxElem):
                maxElem = abs(matrix[k,i])
                maxRow = k
        #Intercambiamos las filas pivotenando el mayor numero en las filas 
        for k in range(i,n+1):
            temp = matrix[maxRow,k]
            matrix[maxRow,k]=matrix[i,k]
            matrix[i,k]=temp

		#Restamos las lineas seleccionadas
        for k in range(i+1,n):
            pivote = -matrix[k,i]/float(matrix[i,i]) #identifico el pivote
            for j in range(i, n+1):
                matrix[k,j] += pivote*matrix[i,j] #Multiplico por el pivote y resto

		#hacemos cero las filas debajo de la columna actual
        for k in range(i+1, n):
            matrix[k,i]=0
    #Resolvemos la matriz triangular superior Ax=b
    A = matrix[:,:m-1]
    b = matrix[:,m-1:m]
    return solve_triang_superior(A,b)
#%%
A_small = crearMatriz("M_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_small = np.matrix(A_small)#transformamos el array en numpy matriz

b_small = crearMatriz("V_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_small = np.matrix(b_small)#transformamos el array en numpy matriz

print("\n---------------Eliminación Gaussiana con pivoteo-------------------\n")
#Para matriz pequeña 
piv=eliminacion_gaussiana_piv(A_small, b_small) #Utilizo eliminacion gaussiana con pivoteo
print("La solución de A=", A_small,"\n y b=",b_small,"usando eliminación gaussiana con pivoteo es:" )
print(piv)
error_piv_small = np.linalg.norm(np.dot(A_small,piv)-np.transpose(b_small))
print("Con un error de: ", error_piv_small)
#%%
#Para matriz grande
A_large = crearMatriz("M_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_large = np.matrix(A_large)#transformamos el array en numpy matriz

b_large = crearMatriz("V_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_large = np.matrix(b_large)#transformamos el array en numpy matriz

piv2=eliminacion_gaussiana_piv(A_large, b_large) #Utilizo eliminacion gaussiana con pivoteo
print("La solución de la matriz 100x100 usando eliminación gaussiana con pivoteo es:" )
print(piv2)
error_piv_large = np.linalg.norm(np.dot(A_large,piv2)-np.transpose(b_large))
print("Con un error de: ", error_piv_large)