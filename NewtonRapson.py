import numpy as np
import math
import matplotlib.pyplot as plt
import random

#funcion
f = lambda x: 204165.5/(330-2*x) + 10400/(x-20)

#primera derivada
fdx = lambda x,dx: (f(x+dx) - f(x-dx))/(2*dx)

#segunda derivada
fddx = lambda x, dx: (f(x+dx)-2*f(x)+f(x-dx))/(dx**2)



def newtonRapson(x0,fdx,fddx,epsilon,max_iter):
    puntoInicial = x0
    print("punto inicial: "+ str(puntoInicial))
    for i in range(0, max_iter):
        fx = f(puntoInicial)
        nuevoPunto = puntoInicial -(fdx(puntoInicial,e)/fddx(puntoInicial,e))

        if(abs(nuevoPunto-puntoInicial)<epsilon):
            return nuevoPunto
            print("terminado num de interaciones:" + str(i+1))

        puntoInicial = nuevoPunto
        print("iteracion no "  + str(i+1))
        print("nuevo punto:"+ str(puntoInicial))

    return  puntoInicial







#temperatura
T = np.arange(40,90.001, 0.001)

#punto inicial
#puntoInicial  = random.uniform(40.0,90.0)
puntoInicial  = 90.0
e = math.exp(1)**-4



fx = f(newtonRapson(puntoInicial,fdx,fddx,e,100))

print(fx)