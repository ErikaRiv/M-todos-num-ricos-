#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 21:02:37 2020

@author: erikarivadeneira
"""
#%%
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt 

#%% INTERPOLACIÓN DE GREGORY-NEWTON
#xi  vector de puntos xi
#y vector de imagenes de xi
#funcion:funcion real
def u_cal(u, n): 
    temp = u; 
    for i in range(1, n): 
        temp = temp * (u - i); 
    return temp; 
  
def factorial(n): 
    f = 1; 
    for i in range(2, n + 1): 
        f *= i; 
    return f; 
def gregory_newton(xi,fi,funcion):
    n = len(xi)
    y = np.zeros((n,n))
    for i in range(n):
        y[i,0]=fi[i]
    for i in range(1, n): 
        for j in range(n - i): 
            y[j,i] = y[j + 1,i - 1] - y[j,i - 1]
      
    # Displaying the forward difference table 
    for i in range(n): 
        print(xi[i], end = "\t")
        for j in range(n - i): 
            print(y[i,j], end = "\t")
        print("")
       
    value= sym.Symbol('x')  
    polinomio = y[0,0] 
    u = (value - xi[0]) / (xi[1] - xi[0]); 
    for i in range(1,n): 
        polinomio = polinomio + (u_cal(u, i) * y[0,i]) / factorial(i)
    simple = polinomio.expand()
    print("\nPolinomio: \nP(x) =", polinomio)
    print("\nForma simplificada: \nP(x) =",simple )
    px = sym.lambdify(value,simple)
    muestras = 50
    a = np.min(xi)
    b = np.max(xi)
    p_xi = np.linspace(a,b,muestras)
    pol_y=px(p_xi)
    plt.plot(p_xi,pol_y,color= "black")
    plt.plot(np.arange(xi[0],xi[n-1],0.00001),funcion,color="royalblue")
    plt.plot(xi,fi, 'o', color="rebeccapurple")
    plt.show()
    return polinomio

#%%--------FUNCIONES A CONSIDERAR-------------
#-----------------FUNCION 1)-------------------------------------------
print("-----FUNCION f(x)=sin(x)-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
punto=-1
x = sym.Symbol('x')
xi = np.arange(-np.pi,np.pi,2)
y = np.sin(xi)
pol1 = gregory_newton(xi,y,np.sin(np.arange(xi[0],xi[len(xi)-1],0.00001)))
px1 = sym.lambdify(x,pol1) 
error = abs(np.sin(-2)-px1(-2))
print("\nError: ", error)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi2 = np.arange(-np.pi,np.pi,1.5)
y2 = np.sin(xi2)
pol2 = gregory_newton(xi2,y2,np.sin(np.arange(xi2[0],xi2[len(xi2)-1],0.00001)))
px2 = sym.lambdify(x,pol2) 
error2 = abs(np.sin(punto)-px2(punto))
print("\nError: ", error2)
#%%
#TERCERA MUESTRA
print("\n---------Tercera muestra ------\n")
xi3 = np.arange(-np.pi,np.pi,1)
y3 = np.sin(xi3)
pol3 = gregory_newton(xi3,y3,np.sin(np.arange(xi3[0],xi3[len(xi3)-1],0.00001)))
px3 = sym.lambdify(x,pol3) 
error3 = abs(np.sin(punto)-px3(punto))
print("\nError: ", error3)
#%%
#-----------------FUNCION 2)-------------------------------------------
print("-----FUNCION f(x)=x2(5x−3)−2x^4+4x−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi4 = np.arange(-2,4.1,3)
y4 = xi4**2*(5*xi4-3)-2*xi4**4+4*xi4-5
xi4_1 = np.arange(-2,4,0.00001)
funcion = ((xi4_1**2)*(5*xi4_1-3)-2*xi4_1**4+4*xi4_1-5)
pol4 = gregory_newton(xi4,y4,funcion)
px4 = sym.lambdify(x,pol4) 
error4 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px4(punto))
print("\nError: ", error4)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi5 = np.arange(-2,4.1,2)
y5 = xi5**2*(5*xi5-3)-2*xi5**4+4*xi5-5
pol5 = gregory_newton(xi5,y5,funcion)
px5 = sym.lambdify(x,pol5) 
error5 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px5(punto))
print("\nError: ", error5)
#%%
#TERCERA MUESTRA
print("\n---------Tercera muestra ------\n")
xi5_1 = np.arange(-2,4.1,1.5)
y5_1 = xi5_1**2*(5*xi5_1-3)-2*xi5_1**4+4*xi5_1-5
pol5_1 = gregory_newton(xi5_1,y5_1,funcion)
px5_1 = sym.lambdify(x,pol5_1) 
error5_1 = abs((punto**2*(5*punto-3)-2*punto**4+4*punto-5)-px5_1(punto))
print("\nError: ", error5_1)
#%%
#-----------------FUNCION 3)-------------------------------------------
print("-----FUNCION f(x)=x^2−5-------\n")
#PRIMERA MUESTRA
print("\n---------Primera muestra ------\n")
xi6 = np.arange(-10,10.1,10)
y6 = xi6**2-5
xi6_1 = np.arange(-10,10,0.00001)
funcion2 = xi6_1**2-5
pol6 = gregory_newton(xi6,y6,funcion2)
px6 = sym.lambdify(x,pol6) 
error6= abs(((2.5)**2-5)-px6(2.5))
print("\nError: ", error6)
#%%
#SEGUNDA MUESTRA
print("\n---------Segunda muestra ------\n")
xi7 = np.arange(-10,10.1,5)
y7 = xi7**2-5
pol7 = gregory_newton(xi7,y7,funcion2)
px7 = sym.lambdify(x,pol7) 
error7= abs((2.5**2-5)-px7(2.5))
print("\nError: ", error7)