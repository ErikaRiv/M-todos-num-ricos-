#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 19:13:02 2020

@author: erikarivadeneira
"""
#%%CARGANDO LIBRERIAS 
import numpy as np
#%% INTEGRACION DE ROMBERG 
#funcion: funcion a considerar
#a: extremo inferior de interavalo de integracion
#b:  extremo superior de interavalo de integracion
#n: iteraciones del método
def romberg_int( funcion, a, b, n ):
    r = np.array( [[0] * (n+1)] * (n+1), float )
    h = b - a
    r[0,0] = 0.5 * h * ( funcion( a ) + funcion( b ) )
    p2 = 1
    for i in range( 1, n + 1 ):
        h = 0.5 * h #medio paso y puntos restantes
        suma = 0
        p2 = 2 * p2
        for k in range( 1, p2, 2 ):
            suma = suma + funcion( a + k * h )
        r[i,0] = 0.5 * r[i-1,0] + suma * h #regla del trapezoide compuesto para el siguiente nivel de subdivisión
        p4 = 1
        for j in range( 1, i + 1 ):
            p4 = 4 * p4
            r[i,j] = r[i,j-1] + ( r[i,j-1] - r[i-1,j-1] ) / ( p4 - 1 )

    return r
#%% FUNCION A CONSIDERAR 
def f(x):
    y = 1/(1+x)
    return y
#%%IMPRESION DE RESULTADOS
print(romberg_int( f, 0, 1, 6 ))