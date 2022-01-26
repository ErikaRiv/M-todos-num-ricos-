#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 20:04:09 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from Matriz_diagonal import crearMatriz
from Solver_triangular_sup import solve_triang_superior
from Solver_triangular_inf import solve_triang_inferior
#%%
#Creamos funcion para el metodo de Crout 
def crout( A_crout ):
    n=len(A_crout)
    #Creamos matrices de ceros para guardar L Y U
    L = np.zeros((n,n)) 
    U = np.zeros((n,n))    
    for j in range(n):
        U[j,j] = 1 #Llenamos de 1 la diagonal de U
        for i in range(j,n): #Llenamos hacia atrás la matriz L 
            Ltemp = A_crout[i,j] #Guardamos datos de operaciones para L en una variable temporal
            for k in range(j):
                Ltemp -= L[i,k]*U[k,j]
            L[i,j] = Ltemp
        for i in range(j+1,n): #Llenamos hacia adelante la matriz U
            Utemp = A_crout[j,i] #Guardamos datos de operaciones para U en una variable temporal
            for k in range(j):
                Utemp -= L[j,k]*U[k,i]
            #En el caso que un elemento de la diagonal de L sea cero la arreglamos con un valor cercano
            if int(L[j,j]) == 0: 
                L[j,j] = 1e-13
            U[j,i] = Utemp/L[j,j]
    return L, U #Retornamos las matrices L U 

#%%
#--------------------Para el uso del método de Crout--------------------------

print("\n------------Descomposición LU por el método de Crout-----------------\n")
#Probamos con la matriz pequeña
A_small = crearMatriz("M_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_small = np.matrix(A_small)#transformamos el array en numpy matriz

b_small = crearMatriz("V_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_small = np.matrix(b_small)#transformamos el array en numpy matriz

Crout_small = crout(A_small)
print("La descomposición LU de\n", A_small, "\n es \n L=", Crout_small[0],"\nU=",Crout_small[1])

L_small = Crout_small[0]#MATRIZ L DE DESCOMPOSICION LU
U_small = Crout_small[1]#MATRIZ U DE DESCOMPOSICION LU

print("\nLU=", np.dot(L_small,U_small))
#Resolvemos el sistema
ytemp = solve_triang_inferior(L_small, b_small)#creamos vector temporal y de Ly=b
sol_LU = solve_triang_superior(U_small, ytemp)#Solucion de Ux=y
print("\n-----------Solución de la descomposición LU---------")
print("La solución del sistema usando factorización LU es: \n",sol_LU)
error_LU=np.linalg.norm(np.dot(A_small,sol_LU)-np.transpose(b_small))#Error de solucion
print("Error:", error_LU,"\n")
#%%
#Probamos con la matriz grande
A_large = crearMatriz("M_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_large = np.matrix(A_large)#transformamos el array en numpy matriz

b_large = crearMatriz("V_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_large = np.matrix(b_large)#transformamos el array en numpy matriz

Crout_large = crout(A_large)

L_large = Crout_large[0]#MATRIZ L DE DESCOMPOSICION LU
U_large = Crout_large[1]#MATRIZ U DE DESCOMPOSICION LU
#Resolvemos el sistema
ytemp2 = solve_triang_inferior(L_large, b_large)#creamos vector temporal y de Ly=b
sol_LU2= solve_triang_superior(U_large, ytemp2)#Solucion de Ux=y
print("La descomposición LU de\n", A_large, "\n es \nL=", Crout_large[0], "\nU=",Crout_large[1])
print("\n-----------Solución de la descomposición LU---------")
print("La solución del sistema usando factorización LU es: \n",sol_LU2)
error_LU2=np.linalg.norm(np.dot(A_large,sol_LU2)-np.transpose(b_large))#Error de solucion
print("Error:", error_LU2,"\n")