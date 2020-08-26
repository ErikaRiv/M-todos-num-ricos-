import numpy as np
import math
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero

#-------------------------CREANDO FIGURA CON COORDENADAS EN PIXELES----------------------------------------------------
#Como la funcion range solo trabaja con enteros, creamos un range para flotantes
def frange(x_inicial, x_final, paso):
    while x_inicial < x_final:
        yield x_inicial
        x_inicial += paso
        
#El dominio de la funcion es (-2pi,2pi)
x= np.array(list(frange(-2*np.pi, 2*np.pi,0.1)))
#Número de puntos de x
n = len(x)
#Consideramos la función y(x)=sin(x)
def f(x):
    y = np.sin(x)
    return y
y = f(x)
#Definimos los límites de x
def coor_screen( x_real, y_real ):
    n = len(x_real)
    a = x_real[0]
    b = x_real[n-1]
    #Encontramos el incremento para x 
    delta_x = (b-a)/n
    #Encontramos el máximo y mínimo de y 
    y_min = 1e30
    y_max = 1e-30
    #print(y)
    for i in range(len(y_real)):
        if y_min > y_real[i]:
            y_min = y_real[i]
        if y_max < y_real[i]:
            y_max = y_real[i]
    #HASTA EL MOMENTO TENEMOS DEFINIDAS NUESTRAS COORDENADAS REALES 
    #Pixeles de mi compu 2560 x 1600
    #Dimensión de la ventana del plot 
    size_h = b - a
    size_v = y_max - y_min
    #Asumiendo que tenemos una pantalla de 640x480 pixeles 
    pixel_h = 640
    pixel_v = 480
    #Cálculo del tamaño de cada pixel
    p_delta_x = (b-a)/(pixel_h-20) #ancho del pixel
    p_delta_y = (y_max-y_min)/(pixel_v-20) #largo del pixel
    px = np.zeros(len(x_real))
    py = np.zeros(len(y_real))
    #Cómputo de las coordenadas en pantalla 
    for i in range(len(x_real)):
        px[i] = (x_real[i]/p_delta_x)+10
    for j in range(len(y)):
        py[j]= ((y_max-y[j])/p_delta_y)+10
    return px, py

z = coor_screen(x,y)

#Para proyectar figura
fig = plt.figure(1)
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

for direccion in ["xzero", "yzero"]:
    ax.axis[direccion].set_axisline_style("-|>")
    ax.axis[direccion].set_visible(True)

for direccion in ["left", "right", "bottom", "top"]:
    ax.axis[direccion].set_visible(False)

ax.plot(z[0], z[1]-250) #Restamos 250 porque interpreta el origen del eje x en la parte superior de la pantalla
#TENER EN CUENTA QUE LOS LABELS EN LOS EJES SON LAS COORDENADAS EN PANTALLA
#∫∫∫∫∫∫∫∫savefig("grafica.png", dpi=72)
plt.show()



#-------------------------CREANDO FIGURA CON EL USO DE MATPLOTLIB---------------------------------------
# Crear una figura de 8x6 puntos de tamaño, 80 puntos por pulgada
plt.figure(figsize=(8, 6), dpi=80)

# Crear una nueva subgráfica en una rejilla de 1x1
plt.subplot(1, 1, 1)

# Graficar la función escogida con una línea continua azul de 1 pixel de grosor
plt.plot(x, y, color="blue", linewidth=1.0, linestyle="-")


# Establecer límites del eje x
plt.xlim(-8.0, 8.0)

# Ticks en x
plt.xticks(np.linspace(-8, 8, 10, endpoint=True))

# Establecer límites del eje y
plt.ylim(-1.0, 1.0)

# Ticks en y
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica1.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

# Mostrar resultado en pantalla, el resultado se muestra con coordenadas reales
plt.show()

#----------------HACIENDO ZOOM EN LAS GRAFICAS----------------------------------------

#ZOOM 1 RAIZ 1
plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x[24:41], y[24:41], color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-4.5, 0)
plt.xticks(np.linspace(-4.5, 0, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica2.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()

#ZOOM 2 RAIZ 2
plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x[55:70], y[55:70], color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-1, 1)
plt.xticks(np.linspace(-1, 1, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica3.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()

#-------------------------METODO DE BISECCION--------------------------------------
#Algoritmo para el método de bisección
def biseccion(x, f, tolerancia, max_iteraciones): 
    counter = 0
    a = x[0]
    n = len(x)
    b = x[n-1]
    if (f(a) * f(b) >= 0): 
        print("No hay raíces en este intervalo\n") 
        return
   
    c = a 
    while ((b-a) >= tolerancia and counter < max_iteraciones): 
  
        #Econtrando el punto medio  
        c = (a+b)/2
   
        # Check if middle point is root 
        if (f(c) == 0.0): 
            break
   
        # Decide the side to repeat the steps 
        if (f(c)*f(a) < 0): 
            b = c 
        else: 
            a = c 
        counter += 1
    return c, f(c), counter

#--------------------------MÉTODO DE NEWTON-RAPHSON--------------------------------------
# Creamos la función para Newton-Raphson
def newton_Raphson(xi,f,f_derivada,diferencia_relativa=1, counter=0):
    #Definimos la presición del método y el número máximo de iteraciones
    prec_goal = 1.e-10
    nmax = 1000 #tolerancia
    while diferencia_relativa > prec_goal and counter < nmax:
        #Calculamos la diferencia relativa
        diferencia_relativa = np.abs(f(xi)/f_derivada(xi))
        #Calculamos el siguiente valor x 
        x1 = xi - f(xi)/f_derivada(xi)
        #Intercambiamos el output con el input para la siguiente iteración
        xi = x1
        #Incrementamos el contador
        counter += 1
    return xi,counter

#-------------------------APLICANDO METODOS A FUNCIONES DADAS--------------------------
#FUNCIONES 
#1) x^2 en (-1,1)
x1= np.array(list(frange(-1, 1,0.01)))
#Consideramos la función y(x)=x^2
def f1(x):
    y = pow(x,2)
    return y
y1 = f1(x1)
#Definimos la derivada de la funcion
def f_derivada1(x):
    return 2*x



#Graficamos x^2 para ver donde se encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x1, y1, color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-2, 2)
plt.xticks(np.linspace(-2, 2, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica4.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()
#APLICANDO EL METODO DE BISECCION
#Viendo que la raíz se encuentra cerca en 0 consideramos un intervalo cercano a este punto
print("Usando el método de Bisección: \n")
biseccion(x1,f1,0.01,1000)

#El output nos dice que no hay raíces ya que no hay un cambio de signo en f(a) y f(b)

#APLICANDO NEWTON-RAPHSON
raiz_newton = newton_Raphson(0.1,f1,f_derivada1)
error_newton1 = abs(0.00000001-f1(raiz_newton[0]))
print("Usando el método de Newton-Raphson: \n")
print("La raíz de f(x)=x^2 es",round(raiz_newton[0],8))
print("Número de iteraciones: ",raiz_newton[1])
print("Error relativo:", error_newton1)

#Graficamos y(x)=sin(x) con sus respectivas raices encontradas encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x1, y1, color="blue", linewidth=1.0, linestyle="-")

plt.plot(raiz_newton[0], f1(raiz_newton[0]), marker="o", color="red")
plt.xlim(-4, 4)
plt.xticks(np.linspace(-4, 4, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica5.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()
#----------------------------------------------------------------------------------------
#2) x^3-2x^2+8x-1 en (-3,3)
x2= np.array(list(frange(-3, 3,0.1)))
#Consideramos la función y(x)=x^3-2x^2+8x-1 
def f2(x):
    y = pow(x,3)-2*pow(x,2)+8*x-1
    return y
y2 = f2(x2)
#Definimos la derivada de la funcion
def f_derivada2(x):
    return 3*pow(x,2)-4*x+8

#Graficamos y(x)=x^3-2x^2+8x-1 para ver donde se encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x2, y2, color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-1, 1)
plt.xticks(np.linspace(-1, 1, 10, endpoint=True))
plt.ylim(-5.0,5.0)
plt.yticks(np.linspace(-5, 5, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica6.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()
#APLICANDO EL METODO DE BISECCION
#Viendo que la raíz se encuentra cerca en 0 consideramos un intervalo cercano a este punto
z2 = biseccion(x2,f2,0.01,1000)
error_relativo_z2 = abs(0.00000001-z2[1])
print("Usando el método de Bisección: \n")
print("La raíz de x^3-2x^2+8x-1 se encuentra en la coordenada: ", z2[0:2])
print("La raíz se encontró en la iteración número", z2[2])
print("El error relativo es: ", error_relativo_z2)
print("\n")

#APLICANDO NEWTON-RAPHSON
raiz_newton2 = newton_Raphson(0.15,f2,f_derivada2)
error_newton2 = abs(0.00000001-f2(raiz_newton2[0]))
print("Usando el método de Newton-Raphson: \n")
print("La raíz de y(x)=x^3-2x^2+8x-1 es", round(raiz_newton2[0],10))
print("Número de iteraciones: ", raiz_newton2[1])
print("El error relativo es: ", error_newton2)
print("\n")
      
#Graficamos y(x)=x^3-2x^2+8x-1 con sus respectivas raices encontradas encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x2, y2, color="blue", linewidth=1.0, linestyle="-")

plt.plot(z2[0], z2[1], marker="o", color="red")
plt.plot(raiz_newton2[0], f2(raiz_newton2[0]), marker="o", color="green")


plt.xlim(-1, 1)
plt.xticks(np.linspace(-1, 1, 10, endpoint=True))
plt.ylim(-5.0, 5.0)
plt.yticks(np.linspace(-5, 5, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica7.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()



#3) f(x)=sin(x) en (-2pi,2pi)
x3= np.array(list(frange(-2*np.pi+0.01, 2*np.pi+0.0001,0.1)))
#Consideramos la función y(x)=sin(x)
def f3(x):
    y = np.sin(x)
    return y
y3 = f3(x3)
#Definimos la derivada de la funcion
def f_derivada3(x):
    return np.cos(x)

#Graficamos y(x)=sin(x) para ver donde se encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x3, y3, color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-7, 7)
plt.xticks(np.linspace(-7, 7, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica8.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()
#APLICANDO EL METODO DE BISECCION
#En todo el intervalo 
z3_t = biseccion(x3,f3,0.01,1000)
error_relativo_z3_t = abs(0.00000001-z3_t[1])
print("\nAplicacando el método de bisección: \n")
print("En todo el intervalo: \n")
print("La raíz de sin(x) se encuentra en la coordenada: ", z3_t[0:2])
print("La raíz se encontró en la iteración número", z3_t[2])
print("El error relativo es: ", error_relativo_z3_t)
print("\n")
#Vemos que solo encuentra una raíz

#Como esta funcion tiene multiples raices, escogemos un intervalo que encierre cada riz para aplicar el metodo de biseccion
z3 = biseccion(x3[30:40],f3,0.01,1000)
error_relativo_z3 = abs(0.00000001-z3[1])
print("\nModificando los intervalos:\n")
print("La primera raíz de sin(x) se encuentra en la coordenada: ", z3[0:2])
print("La raíz se encontró en la iteración número", z3[2])
print("El error relativo es: ", error_relativo_z3)
print("\n")
z3_1 = biseccion(x3[60:75],f3,0.01,1000)
error_relativo_z3_1 = abs(0.00000001-z3_1[1])
print("La segunda raíz de sin(x) se encuentra en la coordenada: ", z3_1[0:2])
print("La raíz se encontró en la iteración número", z3_1[2])
print("El error relativo es: ", error_relativo_z3_1)
print("\n")
z3_2 = biseccion(x3[90:115],f3,0.01,1000)
error_relativo_z3_2 = abs(0.00000001-z3_2[1])
print("La tercera raíz de sin(x) se encuentra en la coordenada: ", z3_2[0:2])
print("La raíz se encontró en la iteración número", z3_2[2])
print("El error relativo es: ", error_relativo_z3_2)
print("\n")

#APLICANDO NEWTON-RAPHSON
#Para la primera raíz
raiz_newton3_1 = newton_Raphson(-3.5,f3,f_derivada3)
error_newton3_1 = abs(0.00000001-f3(raiz_newton3_1[0]))
print("La primera raíz de y(x)=sin(x) es", round(raiz_newton3_1[0],10))
print("Número de iteraciones: ", raiz_newton3_1[1])
print("El error relativo es: ", error_newton3_1)
print("\n")

#Para la segunda raíz
raiz_newton3_2 = newton_Raphson(-0.1,f3,f_derivada3)
error_newton3_2 = abs(0.00000001-f3(raiz_newton3_2[0]))
print("La segunda raíz de y(x)=sin(x) es", round(raiz_newton3_2[0],10))
print("Número de iteraciones: ", raiz_newton3_2[1])
print("El error relativo es: ", error_newton3_2)
print("\n")
#Para la tercera raíz
raiz_newton3_3 = newton_Raphson(3,f3,f_derivada3)
error_newton3_3 = abs(0.00000001-f3(raiz_newton3_3[0]))
print("La tercera raíz de y(x)=sin(x) es", round(raiz_newton3_3[0],10))
print("Número de iteraciones: ", raiz_newton3_3[1])
print("El error relativo es: ", error_newton3_3)
print("\n")

#Graficamos y(x)=sin(x) con sus respectivas raices encontradas encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x3, y3, color="blue", linewidth=1.0, linestyle="-")

plt.plot(z3[0], z3[1], marker="o", color="red")
plt.plot(z3_1[0], z3_1[1], marker="o", color="red")
plt.plot(z3_2[0], z3_2[1], marker="o", color="red")

plt.plot(raiz_newton3_1[0], f3(raiz_newton3_1[0]), marker="o", color="green")
plt.plot(raiz_newton3_2[0], f3(raiz_newton3_2[0]), marker="o", color="green")
plt.plot(raiz_newton3_3[0], f3(raiz_newton3_3[0]), marker="o", color="green")

plt.xlim(-7, 7)
plt.xticks(np.linspace(-7, 7, 10, endpoint=True))
plt.ylim(-1.0, 1.0)
plt.yticks(np.linspace(-1, 1, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica9.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()

#4) ln(x) en [0,2)
x4= np.array(list(frange(0, 2,0.1)))
#Consideramos la función y(x)=x^2
def f4(x):
    y = np.log(x)
    return y
y4 = f4(x4)
def f_derivada4(x):
    return 1/x

#Graficamos x^2 para ver donde se encuentran las raíces

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x4, y4, color="blue", linewidth=1.0, linestyle="-")

plt.xlim(-1, 2)
plt.xticks(np.linspace(-1, 2, 10, endpoint=True))
plt.ylim(-10.0, 3.0)
plt.yticks(np.linspace(-10,3, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica10.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()

#APLICANDO EL METODO DE BISECCION
#Viendo que la raíz se encuentra cerca en 0 consideramos un intervalo cercano a este punto
z4 = biseccion(x4,f4,0.01,1000)
error_relativo_z4 = abs(0.00000001-z4[1])
print("Usando el método de Bisección: \n")
print("La raíz de f(x)=ln(x) se encuentra en la coordenada: ", z4[0:2])
print("La raíz se encontró en la iteración número", z4[2])
print("El error relativo es: ", error_relativo_z4)
print("\n")

#APLICANDO NEWTON-RAPHSON
raiz_newton4 = newton_Raphson(0.8,f4,f_derivada4)
error_newton4 = abs(0.00000001-f4(raiz_newton4[0]))
print("Usando el método de Newton-Raphson: \n")
print("La raíz de y(x)=ln(x) es", round(raiz_newton4[0],10))
print("Número de iteraciones: ", raiz_newton4[1])
print("El error relativo es: ", error_newton4)
print("\n")

#Graficamos y(x)=ln(x) con sus respectivas raices encontradas 

plt.figure(figsize=(8, 6), dpi=80)
plt.subplot(1, 1, 1)
plt.plot(x4, y4, color="blue", linewidth=1.0, linestyle="-")

plt.plot(z4[0], z4[1], marker="o", color="red")
plt.plot(raiz_newton4[0], raiz_newton4[1], marker="o", color="green")

plt.xlim(-1, 2)
plt.xticks(np.linspace(-1, 2, 10, endpoint=True))
plt.ylim(-10.0, 3.0)
plt.yticks(np.linspace(-10, 3, 6, endpoint=True))

# Guardar la figura usando 72 puntos por pulgada
#savefig("grafica11.png", dpi=72)
ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

plt.show()







