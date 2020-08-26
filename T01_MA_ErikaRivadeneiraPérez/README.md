GRAFICADOR Y CEROS DE FUNCIONES
El presente programa está dividido en 3 partes. La primera parte muestra dos métodos de graficación, el primer método transforma coordenadas reales de una función dada en coordenadas de pixeles en pantalla. El segundo método usa netamente la librería matplotlib y platea en coordenadas reales. 

La segunda parte del código muestra 2 métodos para conseguir ceros de funciones. El primer método es de Bisección y el segundo de Newton-Raphson. Por último, se consideran cuatro ejemplos de funciones para las cuales se aplican los dos métodos para raíces de funciones. En cada ejemplo se muestra sus respectivas gráficas para localizar los ceros y ya obteniendo los resultados requeridos se grafican las funciones indicando sus raíces (si es que las hay). 

INSTRUCCIONES DE EJECUCION
El código no requiere datos de entrada para el usuario. Si se desea cambiar los datos para gráficos se debe crear primero una secuencia de puntos para el primer eje
x = np.array(list(frange(-2*np.pi, 2*np.pi,0.1)))
-Después definir la imagen de los puntos:
def f(x):
    y = np.sin(x)
    return y
y = f(x)
-Para graficarlos se ingresan los datos en el siguiente comando
plt.plot(x, y)
Al cual se le pueden agregar otras especificaciones

Las funciones "biseccion" y "newton_Raphson" buscan las raíces de las funciones dadas con los métodos que indican su propio nombre y reciben como parámetros (x,funcion, tolerancia, numero_maximo_iteraciones) y (punto_inicial, funcionarios, derivada_funcion) respectivamente.

Prerequisitos
-Numpy
-Matplotlib
-Pylab

Autora
Erika Rivadeneira Pérez - Matemáticas Aplicadas CIMAT