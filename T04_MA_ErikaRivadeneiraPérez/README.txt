Los presentes códigos "factorizacion.py" y "ecu_calor.py" presentan métodos de factorización (LLT y LDLT) y la solución a la ecuación del calor respectivamente.

-abrir factorizacion.py 
-ejecutar por bloques (no requiere entradas del usuario)

Si el usuarios desea probar los métodos con otra dimensión de la matriz solamente se modifica "n", con n-1 ecuaciones. Por ejemplo: 
	- a=ldlt(CrearMatrizCalor(n)) #FACTORIZACION LDLT
	- b=cholesky(CrearMatrizCalor(n)) #FACTORIZACIÓN CHOLESKY

-abrir ecu_calor.py
-ejecutar por bloques (no requiere entradas del usuario)
	-la función solve_ecuacion_calor( Q, K, phi0, phin, n, L) resulve la ecuación del 	calor, si el usuario desea probar esta función para una matriz con otra dimensión 	solamente debe modificar "n"

NOTA: ecu_calor.py contiene funciones adicionales como nodo_central el cual encuentra la posición del nodo central, la función  error el cual calcula el error ||Ax-b|| y la función plot_temp_nodos la cual grafica los valores de los nodos Φ vs. número de nodos. Es opcional hacer uso de estas funciones 


Autora
Erika Rivadeneira Pérez - Matemáticas Aplicadas CIMAT