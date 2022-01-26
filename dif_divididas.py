#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:45:24 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt 
#%%DIFERENCIAS DIVIDIDAS DE NEWTON
#xi  vector de puntos xi
#y vector de imagenes de xi
#funcion:funcion real
def dif_divididas(xi,y,funcion):
    titulo = ['i','xi','fi']
    n=len(xi)
    pos = np.arange(0,n,1)
    datos = np.concatenate(([pos],[xi],[y]))
    datos = np.transpose(datos)
    dif = np.zeros(shape=(n,n),dtype=float  )
    datos = np.concatenate((datos,dif),axis=1)
    n,m = np.shape(datos)
    diag = n-1
    j= 3
    while(j<m):
        titulo.append("df"+str(j-2))
        i=0
        paso = j-2
        while(i<diag):
            num = datos[i+1,j-1]-datos[i,j-1]
            den =  xi[i+paso]-xi[i]
            datos[i,j]=num/den 
            i+=1
        diag = diag-1
        j+=1
    dif = datos[0,3:]
    n = len(dif)
    x = sym.Symbol('x')
    polinomio = y[0]
    for j in range(1,n):
        fact = dif[j-1]
        term = 1
        for i in range(j):
            term = term *(x-xi[i])
        polinomio+= fact*term
    simple = polinomio.expand()
    print(titulo)
    print(datos)
    print("\nPolinomio: \nP(x) =", polinomio)
    print("\nForma simplificada: \nP(x) =",simple )
    px = sym.lambdify(x,simple)
    muestras = 50
    a = np.min(xi)
    b = np.max(xi)
    p_xi = np.linspace(a,b,muestras)
    pol_y=px(p_xi)
    plt.plot(p_xi,pol_y,color= "black")
    plt.plot(np.arange(xi[0],xi[n-1],0.00001),funcion,color="royalblue")
    plt.plot(xi,y, 'o', color="rebeccapurple")
    plt.show()
    return polinomio
#%%--------FUNCIONES A CONSIDERAR-------------
#-----------------FUNCION 1)-------------------------------------------
print("-----FUNCION f(x)=sin(x)-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
punto=-1
x = sym.Symbol('x')
xi = [-3.14,0,1,3.14]#np.arange(-np.pi,np.pi,2)
y = np.sin(xi)
pol1 = dif_divididas(xi,y,np.sin(np.arange(xi[0],xi[len(xi)-1],0.00001)))
px1 = sym.lambdify(x,pol1) 
error = abs(np.sin(-2)-px1(-2))
print("\nError: ", error)

#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi2 = [-3.14,-2,0,2,3.14]#np.arange(-np.pi,np.pi,1.5)
y2 = np.sin(xi2)
pol2 = dif_divididas(xi2,y2,np.sin(np.arange(xi2[0],xi2[len(xi2)-1],0.00001)))
px2 = sym.lambdify(x,pol2) 
error2 = abs(np.sin(punto)-px2(punto))
print("\nError: ", error2)
#%%
#TERCERA MUESTRA
print("\n---------Tercera muestra ------\n")
xi3 = [-3.14,-2,-1,1,2,3.14]#np.arange(-np.pi,np.pi,1)
y3 = np.sin(xi3)
pol3 = dif_divididas(xi3,y3,np.sin(np.arange(xi3[0],xi3[len(xi3)-1],0.00001)))
px3 = sym.lambdify(x,pol3) 
error3 = abs(np.sin(punto)-px3(punto))
print("\nError: ", error3)
#%%
#-----------------FUNCION 2)-------------------------------------------
print("-----FUNCION f(x)=x2(5x−3)−2x^4+4x−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi4 = np.array([-2,2,4])#np.arange(-2,4.1,2)
y4 = xi4**2*(5*xi4-3)-2*xi4**4+4*xi4-5
xi4_1 = np.arange(-2,4,0.00001)
funcion = ((xi4_1**2)*(5*xi4_1-3)-2*xi4_1**4+4*xi4_1-5)
pol4 = dif_divididas(xi4,y4,funcion)
px4 = sym.lambdify(x,pol4) 
error4 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px4(punto))
print("\nError: ", error4)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi5 = np.array([-2,0,3,4])#np.arange(-2,4.1,1.5)
y5 = xi5**2*(5*xi5-3)-2*xi5**4+4*xi5-5
pol5 = dif_divididas(xi5,y5,funcion)
px5 = sym.lambdify(x,pol5) 
error5 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px5(punto))
print("\nError: ", error5)
#%%
#TERCERA MUESTRA
print("\n---------Tercera muestra ------\n")
xi5_1 = np.array([-2,-0.5,1.5,2.3,4])#np.arange(-2,4.1,1.5)
y5_1 = xi5_1**2*(5*xi5_1-3)-2*xi5_1**4+4*xi5_1-5
pol5_1 = dif_divididas(xi5_1,y5_1,funcion)
px5_1 = sym.lambdify(x,pol5_1) 
error5_1 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px5_1(punto))
print("\nError: ", error5_1)
#%%
#-----------------FUNCION 3)-------------------------------------------
print("-----FUNCION f(x)=x^2−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi6 = np.array([-10,-1,10])#np.arange(-10,10.1,10)]
y6 = xi6**2-5
xi6_1 = np.arange(-10,10,0.00001)
funcion2 = xi6_1**2-5
pol6 = dif_divididas(xi6,y6,funcion2)
px6 = sym.lambdify(x,pol6) 
error6= abs((punto**2-5)-px6(punto))
print("\nError: ", error6)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi7 = np.arange(-10,10.1,5)
y7 = xi7**2-5
pol7 = dif_divididas(xi7,y7,funcion2)
px7 = sym.lambdify(x,pol7) 
error7= abs((punto**2-5)-px7(punto))
print("\nError: ", error7)
#%%
