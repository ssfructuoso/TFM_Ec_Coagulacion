import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from decimal import Decimal

#######################################
#### Obtencion grafica de las distribuciones iniciales para algunos nucleos
########################################

def id(varX, varT, M00):
    M0=2*M00/(2+M00*varT)
    f=M0**2*np.exp(-M0*varX)
    return f

def producto(varX, varT):
    varX=mp.mpf(varX)
    varT=mp.mpf(varT)
    T=mp.mpf(0)
    '''
    def I_bessel(x):
        integrando=lambda theta: np.exp(x*np.cos(theta))*np.cos(theta)
        integral=integrar.quad(integrando, 0, np.pi)[0]
        integral=integral/(np.pi)
        return integral
    '''

    if varT <= 1:
        T=1+varT
    else:
        T=2*varT**(0.5)

    #f=np.exp(-T*varX)*I_bessel(2*varX*varT**(0.5))/(varX**2 * varT**(0.5))
    f=mp.exp(-T*varX)*mp.besseli(1,2*varX*varT**(0.5))/(varX**2 * varT**(0.5))

    return f

def suma(varX, varT):
    f=(2*np.pi)**(-0.5)*np.exp(-varT)*varX**(-1.5)*np.exp(-0.5*varX*np.exp(-2*varT))
    return f


a=0.001
b=1000
t=2
points=np.linspace(a,b,num=100*b)

#ID
fId=[]
for i in points:
    fId.append(id(i,t,1.0))


###Suma
fSuma=[]
for i in points:
    fSuma.append(suma(i,t))


###Producto
fProducto=[]
for i in points:
    fProducto.append(producto(i,t))

###
##Plot
plt.figure(1)
plt.plot(points, fId, label='id')
plt.plot(points, fSuma, label='suma')
plt.plot(points, fProducto, label='producto')
plt.legend()

plt.figure(2)
'''
plt.plot(np.log10(points), np.log10(fId), label='id')
plt.plot(np.log10(points), np.log10(fSuma), label='suma')
plt.plot(np.log10(points), np.log10(fProducto), label='producto')
'''

plt.plot(points, np.log10(fId), label='id')
plt.plot(points, np.log10(fSuma), label='suma')

fProductoList=np.zeros(len(fProducto))
for i in range(0,len(fProducto)):
    fProductoList[i]=float(fProducto[i])
plt.plot(points, np.log10(fProductoList), label='producto')

plt.legend()


plt.show()