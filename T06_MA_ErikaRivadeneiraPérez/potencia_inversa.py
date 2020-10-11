"""
Created on Sun Sep 27 13:31:18 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
from numba import jit
#%%--------FUNCIONES SOLVER DE SISTEMAS DE ECUACIONES--------------------------
def crout( A_crout ):
    n=len(A_crout)
    #Creamos matrices de ceros para guardar L Y U
    L = np.zeros((n,n)) 
    U = np.zeros((n,n))    
    for j in range(n):
        U[j,j] = 1 #Llenamos de 1 la diagonal de U
        for i in range(j,n): #Llenamos hacia atrás la matriz L 
            Ltemp = A_crout[i,j] #Guardamos datos de operaciones para L en una variable temporal
            for k in range(j):
                Ltemp -= L[i,k]*U[k,j]
            L[i,j] = Ltemp
        for i in range(j+1,n): #Llenamos hacia adelante la matriz U
            Utemp = A_crout[j,i] #Guardamos datos de operaciones para U en una variable temporal
            for k in range(j):
                Utemp -= L[j,k]*U[k,i]
            #En el caso que un elemento de la diagonal de L sea cero la arreglamos con un valor cercano
            if int(L[j,j]) == 0: 
                L[j,j] = 1e-13
            U[j,i] = Utemp/L[j,j]
    return L, U #Retornamos las matrices L U 
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
 #Sustitución hacia atrás
@jit(nopython=True)
def solve_triang_superior( matriz, b ):
    n=len(matriz) #Número de elementos de A
    x = np.zeros(n)#Creamos vector de ceros para guardar soluciones
    for i in range(n-1,-1,-1):#recorremos la matriz desde el final hasta el inicio
        sumj = 0 #inicializamos las sumas
        for j in range(i+1,n): 
            sumj += matriz[i,j]*x[j] #Realizamos las sumas de cada iteracion
        x[i] = (b[i]-sumj)/matriz[i,i]
    return x #Retornamos las raices 
#%%-------------------METODO DE LA POTENCIA INVERSA----------------------------
#A : matriz a la que se desea encontrar el autovalor más grande
#tol : Tolerancia aceptada
#max_iter : Número máximo de iteraciones aceptadas   
def Potencia_Inversa(A, max_iter, tol):
    n,m = np.shape(A)
    cont = 0
    lambd = 0
    vo = np.random.rand(n)
    v1 = np.zeros(n)
    error=1
    while(cont<max_iter and error > tol ):
        aux = vo.copy()
        vo = vo/np.linalg.norm(vo)
        L,U= crout(A)
        z = solve_triang_inferior(L,vo)
        v1 = solve_triang_superior(U, z)
        suma = 0
        for i in range(n):
            suma += v1[i]*vo[i]
        vo[:]=v1
        suma = 1/suma
        vo=vo/np.linalg.norm(vo)
        error = abs(lambd-suma)
        lambd = suma
        cont+=1
    vo=aux.copy()
    #print("Autovalor: ", lambd)
    #print("Autovector: ", vo)
    #print("Iteración: ", cont)
    #print("Error: ", error)
    return lambd, vo ,cont,error
#%%----------------METODO DE LA POTENCIA INVERSA CON DEFLACION-----------------
#A : matriz a la que se desea encontrar el autovalor más grande
#tol : Tolerancia aceptada
#max_iter : Número máximo de iteraciones aceptadas   
#r: número de autovalores a considerar 
def Potencia_Inversa_Deflacion(A, max_iter, tol, r ):
    n,m = np.shape(A)
    cont = 0
    lambd = 0
    vo = np.random.rand(n)
    v1 = np.zeros(n)
    valores =[]
    vectores = []
    #Primer valor propio más pequeño
    vo1=Potencia_Inversa(A,max_iter,tol)
    valores.append(vo1[0])
    vectores.append(vo1[1])
    error=1
    for i in range(r-1):
        while(cont<max_iter  ):
            aux = vo.copy()
            for j in vectores:
                temp=0
                temp=np.dot(j,vo)
                vo-=temp*j #Quitando las contribuciones
            vo = vo/np.linalg.norm(vo)
            L,U= crout(A)
            z = solve_triang_inferior(L,vo)
            v1 = solve_triang_superior(U, z)
            suma = 0
            for i in range(n):
                suma += v1[i]*vo[i]
            vo[:]=v1
            suma = 1/suma
            vo=vo/np.linalg.norm(vo)
            error = abs(lambd-suma)
            if( error < tol):
                v1=v1/np.linalg.norm(v1)
                valores.append(lambd)
                vectores.append(vo)
                vo1=[suma,v1,cont].copy()
                break
            else:
                lambd = suma
                cont+=1
        vo=aux.copy()
        #print(vo1)
    
    return valores, vectores ,cont,error



#%% CARGANDO MATRICES
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)#A3x3
A1=np.loadtxt("Eigen_50x50.txt", skiprows=1)#A50x50
A2=np.loadtxt("Eigen_500x500.txt", skiprows=1)#A500x500
A3=np.loadtxt("Eigen_1000x1000.txt", skiprows=1)#A1000X1000
#%% TOLERANCIA E ITERACIONES SEGURIDAS
tol = 1e-7
max_iter=10000
#%%---------POTENCIA INVERSA EN MATRIZ 3x3-------------------------------------
print("------------POTENCIA INVERSA EN MATRIZ 3x3-------------\n")
Vals1=Potencia_Inversa(A, max_iter, tol)
print("Autovalor: ",Vals1[0],"\nAutovector: ",Vals1[1],"\nIteraciones: ", Vals1[2],"\nError: ",Vals1[3])
#%%--------------POTENCIA INVERSA CON DEFLACION EN MATRIZ 3x3------------------
print("\n----------POTENCIA INVERSA CON DEFLACION EN MATRIZ 3x3------------\n")
Vals = Potencia_Inversa_Deflacion(A, max_iter, tol, 2 )
print("Autovalores: ",Vals[0],"\nAutovectores: ",Vals[1],"\nIteraciones: ", Vals[2],"\nError: ",Vals[3])
#%%---------POTENCIA INVERSA EN MATRIZ 50x50------------------------
print("------------POTENCIA INVERSA EN MATRIZ 50x50-------------\n")
Vals2=Potencia_Inversa(A1, max_iter, tol)
print("Autovalor: ",Vals2[0],"\nAutovector: ",Vals2[1],"\nIteraciones: ", Vals2[2],"\nError: ",Vals2[3])
#%%------------POTENCIA INVERSA CON DEFLACION EN MATRIZ 50x50------------------
print("\n----------POTENCIA INVERSA CON DEFLACION EN MATRIZ 50x50------------\n")
Vals3 = Potencia_Inversa_Deflacion(A1, max_iter, tol, 10 )
print("Autovalores: ",Vals3[0],"\nAutovectores: ",Vals3[1],"\nIteraciones: ", Vals3[2],"\nError: ",Vals3[3])
#np.savetxt("Autovals_125x125.txt",R1[0], fmt="%f"	)

