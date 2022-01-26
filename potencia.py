#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 20:27:14 2020

@author: erikarivadeneira
"""
#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%%FUNCIONES ADICIONALES
#Función norma para vectores
def Normalizar(v1, num):
    n = len(v1)
    for i in range(n):
        v1[i]=v1[i]/num
    return v1
#Función producto Matriz vector
def prodMV(A,vo):
    n,m = np.shape(A)
    v1 = np.zeros(n)
    for i in range(n):
        suma = 0
        for j in range(n):
            suma += A[i,j]*vo[j]
        v1[i] = suma
    return v1
#%%--------------------------METODO DE LA POTENCIA-----------------------------
#A : matriz a la que se desea encontrar el autovalor más grande
#tol : Tolerancia aceptada
#max_iter : Número máximo de iteraciones aceptadas
def Potencia(A, tol, max_iter):
    n,m = np.shape(A)
    vo = np.random.rand(n)
    cont = 0
    num = 0
    lambd = 0
    error = 1
    while( cont < max_iter and error > tol):
        v1 = prodMV(A,vo)
        num = np.amax(abs(v1))
        arg_pos = np.where(abs(v1) == num)
        num = v1[arg_pos[0]]
        v1 = Normalizar(v1,num)
        vo[:]=v1
        error= abs(lambd-num)
        lambd = num
        cont += 1
    print("Autovalor más grande: \n", num)
    print("Autovector: \n", v1)
    print("Iteración: ", cont)
    print("Error: ", error)
    return lambd, v1
        
#%%--------------------------POTENCIA CON DEFLACION----------------------------
def Potencia_Deflacion(A, r, max_iter, tol):
    m,n = np.shape(A)
    vo = np.random.rand(m)
    autvals = []
    autvecs = []
    vo = vo / np.linalg.norm(vo)
    z = vo.copy()
    for j in range(r):
        for i in range(max_iter):
            aux = vo.copy()
            for k in range(j):
                a = autvecs[k].dot(vo)
                aux = aux - a*autvecs[k]
            vo = aux.copy()
            vo = vo / np.linalg.norm(vo)
            y = A@vo
            mu = y@vo
            y = y / np.linalg.norm(y)
            error = np.linalg.norm((vo - y)) / np.linalg.norm(vo)
            vo = y.copy()
            if error < tol:
                autvals.append(mu)
                autvecs.append(vo)
               # print("Iteraciones: ",i)
              #  print("Autovalor", mu)
               # print("Error: ", error)
                vo = z.copy()
                break
    print("----------------------------")
    print("Autovalores: ", autvals)
    print("Autovectores: ", autvecs)   
    print("Error: ", error)
    return autvals, autvecs
#%%

#%%----------------------------RESULTADOS--------------------------------------
#Cargando matriz de dimensión n=3
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)
#Cargando matriz de dimensión n=125
A1=np.loadtxt("Eigen_125x125.txt",delimiter= ' ', skiprows=1)
#Sugerencia de tolerancia y máximo número de iteraciones 
tol=0.000001
max_iter=1000
#%%---------------PRUEBA METODO DE POTENCIA PARA A3X3--------------------------
Potencia(A, tol, max_iter)

#%%-------------PRUEBA METODO DE POTENCIA PARA A125X125------------------------
Potencia(A1, tol, max_iter)

#%%----------PRUEBA METODO DE POTENCIA CON DEFLACION PARA A3X3------------------
Potencia_Deflacion(A, 2, max_iter, tol)



#%%---------PRUEBA METODO DE POTENCIA CON DEFLACION  PARA A125X125-------------
#R1 = Potencia_Deflacion(A1, 121, max_iter, tol)
#np.savetxt("Autovals_125x125.txt",R1[0], fmt="%f"	)