import numpy as np
import math

f = lambda x: 204165.5/(330 - 2 * x) + 10400/(x - 20)
fdx = lambda x,dx: (f(x+dx) - f(x-dx))/(2*dx)

def cubic(f,fdx,a,b,iter,epsilon):
    dx = 0.01
    x_min = 0
    # for i in range(1,iter):
    #     alpha = (a+b)/2
    #     fxa = fdx(a,dx)
    #     fdalpha = fdx(alpha,dx)
    #     if(fxa*fdalpha<0):
    #         b=alpha
    #         break
    #     else:
    #         a=alpha

    for i in range(0, iter):

        fda = fdx(a,dx)
        fdb = fdx(b,dx)

        z = 3 * (f(a) - f(b)) / (b - a) + fda + fdb
        w = ((b - a) / np.abs(b - a)) * np.sqrt(z * z - fda * fdb)
        miu = calculateMiu(fda,fdb,w,z)
        print('iter ' + str(i) + '  a:' + str(a) + ' b:' + str(b))

        if miu>0 and miu<1:
            x_min = b - miu * (b - a)
        elif miu<0:
            x_min = b
        else:
            x_min = a

        if(abs(fdx(x_min,dx))<epsilon):
            break

        else:
            if(fda*fdx(x_min,dx)<0):
                b = x_min
            else:
                a = x_min

    print("x* =", str(x_min) + "  f(x*)=" + str(f(x_min)))





def calculateMiu(fda,fdb,w,z):
    miu = (fdb + w - z) / (fdb - fda + 2 * w)
    return  miu





epsilon = math.exp(1) ** -3
cubic(f,fdx, 40, 90,100,epsilon)
