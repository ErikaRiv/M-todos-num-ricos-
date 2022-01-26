#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:35:32 2020

@author: erikarivadeneira
"""

#%%
import numpy as np
#%%

#Creamos la función para resolver una matriz diagonal 
def matriz_diagonal( matriz, b ):
    x = np.zeros(len(matriz))
    #if (prod(np.diag(matriz))!=0): #Verificamos que la matriz no sea singular para empezar a operar
    for i in range(len(matriz)):
        x[i]=b[i]/matriz[i,i] 
   # else:
       # print("La matriz es singular") #Si es singul
    return x
#%%
def crearMatriz(name):
    data = open(name,"r")

    return [np.array(list(map(float,i.strip().split(" ")))) for i in data] #splits a string into a list
def imprimir_sol(solucion,matriz,b):
    print('---------Resolviendo Ax=b------------\n donde:\n ')
    print("A =",matriz)
    print("y b =",b)
    print("\nLa solución a este problema es:\n")
    for i in range(len(solucion)):
        
        print('x'+str(i)+' =',solucion[i]) 
    error = np.linalg.norm(np.dot(matriz,solucion)-np.transpose(b))
    print("Error: ", error)
#%%

#--------Para matriz Diagonal-----------------
print("\n---------------Soluciones para mat. diagonales-------------------\n")

A_diag = crearMatriz("M_DIAG.txt")[1:]#Leemos unicamente los elementos de matriz diagonal usando crearMatriz
A_diag = np.matrix(A_diag) #transformamos el array en numpy matriz

b_diag= crearMatriz("V_DIAG.txt")[1:]#Leemos unicamente los elementos de matriz diagonal usando crearMatriz
b_diag = np.matrix(b_diag)#transformamos el array en numpy matriz

Sol_diag = matriz_diagonal(A_diag,b_diag) #Encontramos la solucion del sistema Ax=b para matrices diagonales

imprimir_sol( Sol_diag, A_diag,b_diag)#Imprimimos solucion