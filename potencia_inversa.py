"""
Created on Sun Sep 27 13:31:18 2020

@author: erikarivadeneira
"""

#%%CARGANDO LIBRERIAS----------------------------------------------------------
import numpy as np
#%%
#Función norma para vectores
def Normalizar(v1, num):
    n = len(v1)
    for i in range(n):
        v1[i]=v1[i]/num
    return v1

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
def Potencia_Inversa(A, tol, max_iter):
    n,m = np.shape(A)
    vo = np.zeros(n)
    vo[0]=1
    q = ((np.transpose(vo) @ A) @ vo)/(np.transpose(vo) @ vo)   
    vo = vo/np.linalg.norm(vo)
    cont = 0
    num = 0
    lambd = 0
    error = 1
    A_2 = A-q*np.eye(n)
    fact = crout(A_2)
    L = fact[0]
    U = fact[1]
    while( cont < max_iter  and error > tol):
        z = solve_triang_inferior(L,vo)
        v1 = solve_triang_superior(U, z)
        num = np.amax(abs(v1)) 
        arg_pos = np.where(abs(v1) == num)
        num = v1[arg_pos[0]]
        v1 = Normalizar(v1,num)
        error = np.linalg.norm(vo-v1)/np.linalg.norm(vo)
        lambd = num
        vo = v1.copy()
        cont +=1
    lambd = (1/lambd)+q
    print("Autovalor más pequeño: \n", lambd)
    print("Autovector: \n", vo)
    print("Iteración: ", cont)
    print("Error: ", error)
    return lambd, vo

#%%
def Jacobi(A,max_iter, tol):
    n,m = np.shape(A)            
    lambds = np.zeros(n)     # initialize eigenvalues
    vectores = np.zeros((n,n)) # initialize eigenvector
    cont = 0
    for i in range(0,n): #Colocando unos en la diagonal de matriz de vectores 
        vectores[i,i] = 1
    while(cont<max_iter):
         suma = 0;   
         for i in range(n): 
             suma += sum(abs(A[i,(i+1):n])) #suma de elementos fuera de la diagonal
         if (suma < tol): # diagonal form reached
              for i in range(0,n):
                  lambds[i] = A[i,i]
              break
         else:
              eps = suma/(n*(n-1)/2.0) #Valor promedio fuera de la diagonal
              for i in range(n-1):       
                   for j in range(i+1,n):  
                       if (abs(A[i,j]) > eps):    #Encontrando el elemento más grande fuera de la diagonal
                           val = A[i,i] - A[j,j]       
                           if (abs(val) < tol): #
                               theta = np.pi/4   #ángulo donde tang se hace cero      
                           else: 
                               theta = 0.5*np.arctan(2.0*A[i,j]/val)  #ángulo theta
                           sine = np.sin(theta)
                           cosine = np.cos(theta)
                           #Empizan las rotaciones 
                           for k in range(i+1,j):
                               aux  = A[i,k]
                               A[i,k] = A[i,k]*cosine + A[k,j]*sine  
                               A[k,j] = A[k,j]*cosine - aux *sine  
                           for k in range(j+1,n):
                               aux  = A[i,k]
                               A[i,k] = A[i,k]*cosine + A[j,k]*sine  
                               A[j,k] = A[j,k]*cosine - aux *sine  
                           for k in range(i):
                               aux  = A[k,i]
                               A[k,i] = A[k,i]*cosine + A[k,j]*sine
                               A[k,j] = A[k,j]*cosine - aux *sine
                           aux = A[i,i]
                           A[i,i] = A[i,i]*cosine*cosine + 2*A[i,j]*cosine*sine +A[j,j]*sine*sine
                           A[j,j] = A[j,j]*cosine*cosine - 2*A[i,j]*cosine*sine +aux *sine*sine  
                           A[i,j] = 0                                           
                           for k in range(n):
                                aux  = vectores[k,j]
                                vectores[k,j] = vectores[k,j]*cosine - vectores[k,i]*sine  
                                vectores[k,i] = vectores[k,i]*cosine + aux *sine                                 
         cont+=1
         
    print("Autovalores:\n")
    print(lambds)
    print("\nAutovectores:\n")
    print(vectores)
    print("Iteraciones: ", cont)
    print("Error: ",suma)

    return lambds,vectores,cont

#%%
A=np.loadtxt("Eigen_3x3.txt",delimiter= ' ', skiprows=1)
#%%
A1=np.loadtxt("Eigen_50x50.txt",delimiter= ' ', skiprows=1)
