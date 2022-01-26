#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:50:11 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
from Matriz_diagonal import crearMatriz, imprimir_sol
#%%
#Creamos funcion para resolver matrices triangulares inferiores
def solve_triang_inferior( matriz, b ):
    if( np.allclose(matriz, np.tril(matriz))==True): #Checkando que la matriz sea triangular inferior
        x = np.zeros(len(matriz)) #Creando vector para guardar soluciones
        for i in range(len(matriz)):   #Recorriendo las filas 
            sumj=0 #Inicializando las sumas
            for j in range(i): #Recorriendo columnas
                sumj += matriz[i,j]*x[j] #Realizando la suma de cada iteracion
            x[i]=(b[i]-sumj)/matriz[i,i] 
        return x #Retornando resultado
    else:
        print("La matriz no es triangular inferior")
    
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
#--------Para matriz Triangular Inferior-------------
print("\n---------------Soluciones para mat. triangular inferior-------------------\n")

A_Tinf = crearMatriz("M_TINF.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
A_Tinf = np.matrix(A_Tinf) #transformamos el array en numpy matriz

b_Tinf= crearMatriz("V_TINF.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
b_Tinf = np.matrix(b_Tinf)#transformamos el array en numpy matriz

Sol_Tinf = solve_triang_inferior(A_Tinf,b_Tinf) #Encontramos la solucion del sistema Ax=b para matrices triangular inferior

imprimir_sol( Sol_Tinf, A_Tinf,b_Tinf)#Imprimimos solucion