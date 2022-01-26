#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:10:54 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt 
#%%---------------INTERPOLACIÓN DE LAGRANGE------------------------------------ 
#Esta función tiene incorporada la gráfica de los datos junto al polinomio interpolador
#DATOS DE ENTRADA:
#xi: Vector de datos x
#y: imagenes de xi
#funcion: f(x)
def lagrange(xi,y,funcion):
    n = len(xi)
    x = sym.Symbol('x')
    polinomio = 0 
    #error=0
    for i in range(n):
        num = 1
        den = 1
        for j in range(n):
            if(i!=j):
                num = num * (x-xi[j])
                den = den * (xi[i]-xi[j])
            lagrange = (num/den) * y[i]
        polinomio += lagrange
    polinomio_simple = sym.expand(polinomio)
    px = sym.lambdify(x,polinomio) 
    #error = funcion(punto)-px(punto)
    print("La forma del polinomio de interpolación es:\n\ny =", polinomio )
    print("\nLa forma simplificada del polinomio es:\n\n y =", polinomio_simple)
    print("\nMuestras: ", n)
    #print("\nError: ", error)
    a= np.min(xi)
    b=np.max(xi)
    muestras=20
    pol_xi = np.linspace(a,b,muestras)
    pol_y = px(pol_xi) 
    plt.plot(pol_xi,pol_y,color= "black")
    plt.plot(np.arange(xi[0],xi[n-1],0.00001),funcion,color="royalblue")
    plt.plot(xi,y, 'o', color="rebeccapurple")
    plt.show()
    return polinomio

#def error_interpolacion():
#%%--------FUNCIONES A CONSIDERAR-------------
#-----------------FUNCION 1)-------------------------------------------
print("-----FUNCION f(x)=sin(x)-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
punto=-1
x = sym.Symbol('x')
xi = np.arange(-np.pi,np.pi,1.5)
y = np.sin(xi)
pol1 = lagrange(xi,y,np.sin(np.arange(xi[0],xi[len(xi)-1],0.00001)))
px1 = sym.lambdify(x,pol1) 
error = abs(np.sin(punto)-px1(punto))
print("\nError: ", error)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi2 = np.arange(-np.pi,np.pi,1)
y2 = np.sin(xi2)
pol2 = lagrange(xi2,y2,np.sin(np.arange(xi2[0],xi2[len(xi2)-1],0.00001)))
px2 = sym.lambdify(x,pol2) 
error2 = abs(np.sin(punto)-px2(punto))
print("\nError: ", error2)
#%%
#TERCERA MUESTRA
print("\n---------Tercera muestra ------\n")
xi3 = np.arange(-np.pi,np.pi,0.5)
y3 = np.sin(xi3)
pol3 = lagrange(xi3,y3,np.sin(np.arange(xi3[0],xi3[len(xi3)-1],0.00001)))
px3 = sym.lambdify(x,pol3) 
error = abs(np.sin(punto)-px3(punto))
print("\nError: ", error)
#%%
#-----------------FUNCION 2)-------------------------------------------
print("-----FUNCION f(x)=x2(5x−3)−2x^4+4x−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi4 = np.arange(-2,4.1,2)
y4 = xi4**2*(5*xi4-3)-2*xi4**4+4*xi4-5
xi4_1 = np.arange(-2,4,0.00001)
funcion = ((xi4_1**2)*(5*xi4_1-3)-2*xi4_1**4+4*xi4_1-5)
pol4 = lagrange(xi4,y4,funcion)
px4 = sym.lambdify(x,pol4) 
error4 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px4(punto))
print("\nError: ", error4)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi5 = np.arange(-2,4.1,1.5)
y5 = xi5**2*(5*xi5-3)-2*xi5**4+4*xi5-5
pol5 = lagrange(xi5,y5,funcion)
px5 = sym.lambdify(x,pol5) 
error5 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px5(punto))
print("\nError: ", error5)
#%%
#-----------------FUNCION 3)-------------------------------------------
print("-----FUNCION f(x)=x^2−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi6 = np.arange(-10,10.1,10)
y6 = xi6**2-5
xi6_1 = np.arange(-10,10,0.00001)
funcion2 = xi6_1**2-5
pol6 = lagrange(xi6,y6,funcion2)
px6 = sym.lambdify(x,pol6) 
error6= abs((punto**2-5)-px6(punto))
print("\nError: ", error6)
#%%
#SEGUNDA MUESTRA
print("\n---------Primera muestra ------\n")
xi7 = np.arange(-10,10.1,5)
y7 = xi7**2-5
pol7 = lagrange(xi7,y7,funcion2)
px7 = sym.lambdify(x,pol7) 
error7= abs((punto**2-5)-px7(punto))
print("\nError: ", error7)
#%%
