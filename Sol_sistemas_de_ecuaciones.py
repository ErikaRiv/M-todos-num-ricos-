# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
import numpy as np
#%%

#Creamos la función para resolver una matriz diagonal 
def matriz_diagonal( matriz, b ):
    x = np.zeros(len(matriz))
    #if (prod(np.diag(matriz))!=0): #Verificamos que la matriz no sea singular para empezar a operar
    for i in range(len(matriz)):
        x[i]=b[i]/matriz[i,i] 
   # else:
       # print("La matriz es singular") #Si es singul
    return x


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
#%%

#Creamos función para eliminacion gaussiana sin pivoteo 
def eliminacion_gaussiana(A, b):
    matrix = np.insert(A,A.shape[1],np.transpose(b),1) #Concateno A y b para operar simultaneamente
    np.asarray(matrix) #aseguro el arreglo
    matrix = matrix.astype(float) #hago flotantes a los datos del arreglo
    n,m = matrix.shape
    #empezamos la eliminación
    for i in range(0,n):#row
        for j in range(i+1,n):
            if matrix[j,i] != 0.0:
                #operamos
                multiplier = matrix[j,i]/matrix[i,i]
                matrix[j,i:m]=matrix[j,i:m] - multiplier*matrix[i,i:m] 
    #Resolvemos Ax=b
    A = matrix[:,:n]
    b = matrix[:,n:n+1]
    return solve_triang_superior(A,b)


#%%
#Creamos función para eliminación gaussiana con pivoteo
def eliminacion_gaussiana_piv(A, b):
    matrix = np.insert(A,A.shape[1],np.transpose(b),1)
    n,m = np.shape(matrix) #numero de filas

	# Create an upper triangular matrix
    for i in range(0,n): # for i in [0,1,2,..,n]
		#Encontramos el valor más grande en las columnas 
        maxElem = abs(matrix[i,i])
        maxRow = i
        for k in range(i+1, n):
            if( abs(matrix[k,i])>maxElem):
                maxElem = abs(matrix[k,i])
                maxRow = k
        #Intercambiamos las filas pivotenando el mayor numero en las filas 
        for k in range(i,n+1):
            temp = matrix[maxRow,k]
            matrix[maxRow,k]=matrix[i,k]
            matrix[i,k]=temp

		#Restamos las lineas seleccionadas
        for k in range(i+1,n):
            pivote = -matrix[k,i]/float(matrix[i,i]) #identifico el pivote
            for j in range(i, n+1):
                matrix[k,j] += pivote*matrix[i,j] #Multiplico por el pivote y resto

		#hacemos cero las filas debajo de la columna actual
        for k in range(i+1, n):
            matrix[k,i]=0
    #Resolvemos la matriz triangular superior Ax=b
    A = matrix[:,:m-1]
    b = matrix[:,m-1:m]
    return solve_triang_superior(A,b)
#%%
#----------DESCOMPOSICION LU CON EL METODO DE CROUT---------------------------
#Creamos funcion para el metodo de Crout 
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
#%%
#Creamos funcion para imprimir la solucion de la matriz 

def imprimir_sol(solucion,matriz,b):
    print('---------Resolviendo Ax=b------------\n donde:\n ')
    print("A =",matriz)
    print("y b =",b)
    print("\nLa solución a este problema es:\n")
    for i in range(len(solucion)):
        
        print('x'+str(i)+' =',solucion[i]) 
    error = np.linalg.norm(np.dot(matriz,solucion)-np.transpose(b))
    print("Error: ", error)
#Esta funcion solamente nos sive para imprimir matrices diagonales, triangulares y eliminacion gaussiana sin pivoteo 

#FUNCION PARA LEER ARCHIVOS Y RETORNAR UN NUMPY ARRAY
def crearMatriz(name):
    data = open(name,"r")

    return [np.array(list(map(float,i.strip().split(" ")))) for i in data] #splits a string into a list

#Checkar si las matrices son Diagonal, triangular sup o triangular inf
def tipo_matriz(mat):
    inf = np.allclose(mat, np.tril(mat)) # Para ver si es triangular inferior
    sup = np.allclose(mat, np.triu(mat)) # Para ver si es triangular superior
    diag = np.allclose(mat, np.diag(np.diag(mat))) # Para ver si la matriz es diagonal
    if inf == True:
        print("La matriz es triangular inferior")
    elif sup == True:
        print("La matriz es triangular superior")
    elif diag == True:
        print("La matriz es diagonal")
    else:
        print("La matriz no es diagonal, triangular superior o inferior ")
        
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

#--------Para matriz Diagonal-----------------
print("\n---------------Soluciones para mat. diagonales-------------------\n")

A_diag = crearMatriz("M_DIAG.txt")[1:]#Leemos unicamente los elementos de matriz diagonal usando crearMatriz
A_diag = np.matrix(A_diag) #transformamos el array en numpy matriz

b_diag= crearMatriz("V_DIAG.txt")[1:]#Leemos unicamente los elementos de matriz diagonal usando crearMatriz
b_diag = np.matrix(b_diag)#transformamos el array en numpy matriz

Sol_diag = matriz_diagonal(A_diag,b_diag) #Encontramos la solucion del sistema Ax=b para matrices diagonales

imprimir_sol( Sol_diag, A_diag,b_diag)#Imprimimos solucion
#%%
#--------Para matriz Triangular Inferior-------------
print("\n---------------Soluciones para mat. triangular inferior-------------------\n")

A_Tinf = crearMatriz("M_TINF.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
A_Tinf = np.matrix(A_Tinf) #transformamos el array en numpy matriz

b_Tinf= crearMatriz("V_TINF.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
b_Tinf = np.matrix(b_Tinf)#transformamos el array en numpy matriz

Sol_Tinf = solve_triang_inferior(A_Tinf,b_Tinf) #Encontramos la solucion del sistema Ax=b para matrices triangular inferior

imprimir_sol( Sol_Tinf, A_Tinf,b_Tinf)#Imprimimos solucion
#%%
#--------Para matriz Triangular Superior-------------
print("\n---------------Soluciones para mat. triangular superior-------------------\n")

A_Tsup = crearMatriz("M_TSUP.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
A_Tsup = np.matrix(A_Tsup) #transformamos el array en numpy matriz

b_Tsup= crearMatriz("V_TSUP.txt")[1:]#Leemos unicamente los elementos de matriz triangular usando crearMatriz
b_Tsup= np.matrix(b_Tsup)#transformamos el array en numpy matriz

Sol_Tsup = solve_triang_superior(A_Tsup,b_Tsup) #Encontramos la solucion del sistema Ax=b para matrices triangular superior
imprimir_sol( Sol_Tsup, A_Tsup,b_Tsup)#Imprimimos solucion
#%%
#--------------Para eliminación Gaussiana---------------------
#----------sin pivoteo-------------------------------------------
print("\n---------------Eliminación Gaussiana sin pivoteo-------------------\n")
#Para matriz pequeña
A_small = crearMatriz("M_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_small = np.matrix(A_small)#transformamos el array en numpy matriz

b_small = crearMatriz("V_SMALL.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_small = np.matrix(b_small)#transformamos el array en numpy matriz

Sol_small = eliminacion_gaussiana(A_small,b_small) #Encontramos la solucion del sistema Ax=b para la matriz
imprimir_sol( Sol_small, A_small,b_small)#Imprimimos solucion
#%%
#Para matriz grande
A_large = crearMatriz("M_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
A_large = np.matrix(A_large)#transformamos el array en numpy matriz

b_large = crearMatriz("V_LARGE.txt")[1:]#Leemos unicamente los elementos de matriz  usando crearMatriz
b_large = np.matrix(b_large)#transformamos el array en numpy matriz

Sol_large = eliminacion_gaussiana(A_large,b_large) #Encontramos la solucion del sistema Ax=b para la matriz
imprimir_sol( Sol_large, A_large,b_large)#Imprimimos solucion
#%%
print("\n---------------Eliminación Gaussiana con pivoteo-------------------\n")
#Para matriz pequeña 
piv=eliminacion_gaussiana_piv(A_small, b_small) #Utilizo eliminacion gaussiana con pivoteo
print("La solución de A=", A_small,"\n y b=",b_small,"usando eliminación gaussiana con pivoteo es:" )
print(piv)
error_piv_small = np.linalg.norm(np.dot(A_small,piv)-np.transpose(b_small))
print("Con un error de: ", error_piv_small)
#%%
#Para matriz grande
piv2=eliminacion_gaussiana_piv(A_large, b_large) #Utilizo eliminacion gaussiana con pivoteo
print("La solución de la matriz 100x100 usando eliminación gaussiana con pivoteo es:" )
print(piv2)
error_piv_large = np.linalg.norm(np.dot(A_large,piv2)-np.transpose(b_large))
print("Con un error de: ", error_piv_large)
#%%
#--------------------Para el uso del método de Crout--------------------------

print("\n------------Descomposición LU por el método de Crout-----------------\n")
#Probamos con la matriz pequeña
Crout_small = crout(A_small)
print("La descomposición LU de\n", A_small, "\n es \n L=", Crout_small[0],"\nU=",Crout_small[1])

L_small = Crout_small[0]#MATRIZ L DE DESCOMPOSICION LU
U_small = Crout_small[1]#MATRIZ U DE DESCOMPOSICION LU

print("\nLU=", np.dot(L_small,U_small))
#Resolvemos el sistema
ytemp = solve_triang_inferior(L_small, b_small)#creamos vector temporal y de Ly=b
sol_LU = solve_triang_superior(U_small, ytemp)#Solucion de Ux=y
print("\n-----------Solución de la descomposición LU---------")
print("La solución del sistema usando factorización LU es: \n",sol_LU)
error_LU=np.linalg.norm(np.dot(A_small,sol_LU)-np.transpose(b_small))#Error de solucion
print("Error:", error_LU,"\n")
#%%
#Probamos con la matriz grande
Crout_large = crout(A_large)

L_large = Crout_large[0]#MATRIZ L DE DESCOMPOSICION LU
U_large = Crout_large[1]#MATRIZ U DE DESCOMPOSICION LU
#Resolvemos el sistema
ytemp2 = solve_triang_inferior(L_large, b_large)#creamos vector temporal y de Ly=b
sol_LU2= solve_triang_superior(U_large, ytemp2)#Solucion de Ux=y
print("La descomposición LU de\n", A_large, "\n es \nL=", Crout_large[0], "\nU=",Crout_large[1])
print("\n-----------Solución de la descomposición LU---------")
print("La solución del sistema usando factorización LU es: \n",sol_LU2)
error_LU2=np.linalg.norm(np.dot(A_large,sol_LU2)-np.transpose(b_large))#Error de solucion
print("Error:", error_LU2,"\n")