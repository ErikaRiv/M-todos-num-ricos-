SOLUCIONES Y FACTORIZACIÓN DE SISTEMAS DE ECUACIONES
El presente código en su primera parte contiene funciones para resolver y factorizar distintos istemas de ecuaciones para casos de 
-matrices diagonales, matriz_diagonal( A, b ), output: solucion
-matrices triangulares inferiores (sustitución hacia adelante), solve_triang_inferior( A, b ), output: solucion
-matrices triangulares superiores (sustitución hacia atrás), solve_triang_superior( A, b ),output: solucion
Y métodos como
-eliminación gaussiana, eliminacion_gaussiana(A, b), output: solucion
-función para realizar pivote, pivoteo(A,b), output: solucion
-eliminación gaussiana con pivotes, eliminacion_gaussiana_piv(A, b), output: solucion
-método para factorización LU de Crout, crout( A ), output: L, U
Adicionalmente se encuentran funciones para comprobar el tipo de matriz, imprimir soluciones y leer archivos de texto con una matriz de retorno.

En la segunda parte del programa se implementan las funciones descritas con matrices que fueron leídas como archivos de texto txt. 

El código no solicita ningún input al usuario. Si se desea probar alguna función con otra matriz en archivo txt se debe utilizar la función crearMatriz("file.txt")

Prerequisitos
-Numpy
-Pylab

Autora
Erika Rivadeneira Pérez - Matemáticas Aplicadas CIMAT