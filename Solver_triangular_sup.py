#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:57:03 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from Matriz_diagonal import crearMatriz, imprimir_sol
#%%
#Creamos funcion para resolver matrices triangulares superiores
def solve_triang_superior( matriz, b ):
    n=len(matriz) #Número de elementos de A
    x = np.zeros(n)#Creamos vector de ceros para guardar soluciones
    for i in range(n-1,-1,-1):#recorremos la matriz desde el final hasta el inicio
        sumj = 0 #inicializamos las sumas
        for j in range(i+1,n): 
            sumj += matriz[i,j]*x[j] #Realizamos las sumas de cada iteracion
        x[i] = (b[i]-sumj)/matriz[i,i]
    return x #Retornamos las raices 
#%%
#------------------IMPREMENTACION CON MATRICES.TXT-----------------------------

#¡¡¡¡¡¡¡¡¡¡IMPORTANTE!!!!!!!!!!!!!!


#PARA QUE LA FUNCION LEA LAS MATRICES EN TEXTO LOS NUMEROS CADA ENTRADA DE LA MATRIZ 
#DEBE CONTENER SOLAMENTE UN ESPACIO CON RESPECTO A OTRA ES DECIR SOLO LEERA MATRICES DE ESTE TIPO

#1.O 0.0 0.0 
#0.0 2.0 0.0
#O.O O.O 10.O 

#Y NO LEERÁ

#1.O 0.0  0.0 
#0.0 2.0  0.0
#O.O O.O 10.O


#!!!!!!!!FIN DEL COMUNICADO!!!!!!!!!!!!!!!!!!! 
#%%
#--------Para matriz Triangular Superior-------------
print("\n---------------Soluciones para mat. triangular superior-------------------\n")

A_Tsup = crearMatriz("M_TSUP.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
A_Tsup = np.matrix(A_Tsup) #transformamos el array en numpy matriz

b_Tsup= crearMatriz("V_TSUP.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
b_Tsup= np.matrix(b_Tsup)#transformamos el array en numpy matriz

Sol_Tsup = solve_triang_superior(A_Tsup,b_Tsup) #Encontramos la solucion del sistema Ax=b para matrices triangular superior
imprimir_sol( Sol_Tsup, A_Tsup,b_Tsup)#Imprimimos solucion
