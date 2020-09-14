#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 16:16:31 2020

@author: erikarivadeneira
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from factorizacion import cholesky
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
#------------------SOLVER ECUACION DEL CALOR-----------------------------------
def solve_ecuacion_calor( Q, K, phi0, phin, n, L):
    #Creo matriz tridiagonal
    deltax = L/n
    n=n-1 #Número de ecuaciones
    b = np.ones(n) * (Q*deltax*deltax)/K
    b[0] = b[0]+phi0 #Condiciones Dirichlet
    b[n-1] = b[n-1]+phin
    dl = -1 * np.ones(n)
    du = -1 * np.ones(n)
    d0 = 2 * np.ones(n)
    d = np.vstack((dl, d0, du))
    A = sparse.spdiags(d, (-1, 0, 1), n , n )
    A = A.todense()
    Matrix = np.copy(A)
    L=cholesky(A)
    l = L[0]
    lt = L[1]
    y = solve_triang_inferior(l,b)
    x = solve_triang_superior(lt,y)
    sol = np.copy(x)
    x = np.transpose(np.matrix(x))
    b=np.matrix(b)
    err = error(Matrix,x,b)
    print("n = ",n+1,"\nSolución: ",np.transpose(x),"\nError:",err)
    return sol

#%%
#----------------------FUNCIONES ADICIONALES----------------------------------
#Encuentra la posición del nodo central
def nodo_central(nodos):
    central = int((len(nodos)-1)/2)
    return central 
#CALCULA EL ERROR ||Ax-b||
def error(A, x, b):
    error = np.linalg.norm(np.dot(A,x)-np.transpose(b))
    return error
#Plotea Valores de los nodos Φ vs. número de nodo
def plot_temp_nodos(sol_eq_calor):
    t100 = np.append(10,sol_eq_calor)
    t100 = np.append(sol_eq_calor,20)
    plt.plot(np.array([i for i in range(1,len(t100)+1)]),t100,'c',marker='o',mfc='red')
    plt.xlabel('Nodo')
    plt.ylabel('Temperatura')
    #title("Variación de temperatura en",len(sol_eq_calor)+2 ,"nodos")
#%%
#------------SOLUCION ECUACION DEL CALOR CONSIDERANDO N=4----------------------

n4 = solve_ecuacion_calor( 3, 5, 10, 20, 4, 1)

#%%---------SOLUCION ECUACION DEL CALOR CONSIDERANDO N=100--------------------

n100 = solve_ecuacion_calor( 3, 5, 10, 20, 100, 1)

#------------------------------------------------------------------------------
#------------HASTA AQUI LO REQUERIDO EN LA TAREA-------------------------------
#------------------------------------------------------------------------------

#SOLAMENTE SI EL USUARIO DESEA SE PUEDE EJECUTAR ESTA SECCION
#%%
#------------------------CONVERGENCIA------------------------------------------
n10 = solve_ecuacion_calor( 3, 5, 10, 20, 10, 1)
n30 = solve_ecuacion_calor( 3, 5, 10, 20, 30, 1)
n50 = solve_ecuacion_calor( 3, 5, 10, 20, 50, 1)
n70 = solve_ecuacion_calor( 3, 5, 10, 20, 70, 1)
x1=n4[nodo_central(n4)]
x2=n10[nodo_central(n10)]
x3=n30[nodo_central(n30)]
x4=n50[nodo_central(n50)]
x5=n70[nodo_central(n70)]
x6=n100[ nodo_central(n100)]
x=np.array((x1,x2,x3,x4,x5,x6))
y=np.array((4,10,30,50,70,100))
plt.plot(y,np.transpose(x),'c',marker='o',mfc='red')
plt.xlabel('Nodos')
plt.ylabel('Temperatura')
plt.ylim(15.07, 15.08)  
plt.show()

#%%PLOT NODOS vs TEMPERATURA
plot_temp_nodos(n10)
#%% PLOT DE LA BARRA
plt.imshow(np.matrix(n100),cmap = 'hot', interpolation = 'nearest')
plt.xlabel('Nodos')
plt.yticks([])
plt.ylabel('Temperatura')
