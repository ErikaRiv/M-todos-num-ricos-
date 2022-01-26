#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:39:06 2020

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
#A : matriz a la que se desea encontrar el autovalor más grande
#r: número de autovalores a encontrar 
#tol : Tolerancia aceptada
#max_iter : Número máximo de iteraciones aceptadas
def Potencia_Deflacion(A, r, tol, max_iter):
    n,m = np.shape(A)
    vo = np.random.rand(n)
    w = [] #Para almacenamiento de autovalores
    vp = [] #Para almacenamiento de autovectores
    vo = vo / np.linalg.norm(vo)
    aux = vo.copy()
    for j in range(r):
        for i in range(max_iter):
            aux = vo.copy()
            for k in range(j):
                a = vp[k].dot(vo)
                aux = aux - a*vp[k]
            vo = aux.copy()
            vo = vo / np.linalg.norm(vo)
            y = A@vo
            lambd = y@vo
            y = y / np.linalg.norm(y)
            error = np.linalg.norm((vo - y)) / np.linalg.norm(vo)
            vo = y.copy()
            if error < tol:
                w.append(lambd)
                vp.append(vo)
                print("Iteración: ", i+1)
                print("Autovalor: ", lambd)
                vo = aux.copy()
                break
    print("Autovalores: ", w)
    print("Autovectores: ", vp)
    return w, vp

def Potencia_Deflacion2(A, r, tol, max_iter):
    valores = []
    vectores = []
    n,m = np.shape(A)
    vo = np.random.rand(n)
    lambd = 0
    error = 1
    cont = 0
    for i in range(r):       
        while(error < tol and cont <max_iter):
            resta = vo.copy() 
            for j in range(i):
                a = vectores[i].dot(vo)
                resta = resta - a * vectores[j]
            vo = resta.copy()
            vo = vo / np.linalg.norm(vo)
            v1= prodMV(A, vo)
            lambd = v1@vo
            num = 0
            for k in range(m):
                num += v1[k]*vo[k]
                vo[k]=v1[k]
            vo = vo / np.linalg.norm(vo)
            error = abs(num-lambd)
            v1 = v1/np.linalg.norm(v1)
            valores.append(lambd)
            vectores.append(v1)
            lambd = num
            cont+=1
    return valores, vectores 

def Potencia_2(A,r, tol, max_iter):
    n,m = np.shape(A)
    valores = np.ones(r)
    vectores = []
    vo = np.random.rand(n)
    cont = 0
    num = 0
    lambd = 0
    error = 1
    resta = vo.copy()
    for i in range(r):
        while( cont < max_iter ):
            resta = vo.copy() 
            for j in range(i):
                a = vectores[i].dot(vo)
                resta = resta - a * vectores[j]
            vo = resta.copy()
            v1 = prodMV(A,vo)
            num = np.amax(abs(v1))
            arg_pos = np.where(abs(v1) == num)
            num = v1[arg_pos[0]]
            v1 = Normalizar(v1,num)
            vo[:]=v1
            error= abs(lambd-num)
            if(error < tol):
                valores[i]=lambd
                vectores.append(v1)
                break 
            lambd = num
            cont += 1
       # print("Autovalor más grande: \n", num)
        #print("Autovector: \n", v1)
        #print("Iteración: ", cont)
        #print("Error: ", error)
    return valores, vectores
            
#%%
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)
#Cargando matriz de dimensión n=125
A1=np.loadtxt("Eigen_125x125.txt",delimiter= ' ', skiprows=1)
#Sugerencia de tolerancia y máximo número de iteraciones 
tol=0.000001
max_iter=1000
#%%---------------PRUEBA METODO DE POTENCIA PARA A3X3--------------------------
print("--------Método de Potencia---------")
Potencia(A, tol, max_iter)

#%%-------------PRUEBA METODO DE POTENCIA PARA A125X125------------------------
print("--------Método de Potencia---------")
Potencia(A1, tol, max_iter)

#%%--------PRUEBA METODO DE POTENCIA CON DEFLACION PARA A3X3-------------------
print("--------Método de Potencia con Deflación---------")
Potencia_Deflacion(A, 2,tol, max_iter)

#%%-------PRUEBA METODO DE POTENCIA CON DEFLACION PARA A125X125----------------
print("--------Método de Potencia con Deflación---------")
Potencia_Deflacion(A, 121,tol, max_iter)