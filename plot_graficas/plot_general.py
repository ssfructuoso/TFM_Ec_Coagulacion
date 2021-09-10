import numpy as np
import matplotlib.pyplot as plt
import pandas
import tailer as tl
import io
from decimal import Decimal

###################################
### Representacion grafica de las aproximaciones obtenidas en los casos de Oceanografia
### Figuras del capitulo 3 de la memoria del trabajo
###################################


datosTxt='ficheros_simulaciones/SECCIONAL/RK3/trapecios_compuesto/seccional_RK3_ID_trapecios_compuesto_100.txt'


dfHeader = pandas.read_csv(datosTxt, sep=';',delimiter=';',nrows=1,names=["0","1","2","3","4","5","6"])

Nx=dfHeader.iloc[0,0]
iniX12=dfHeader.iloc[0,1]
R=dfHeader.iloc[0,2]
t0=dfHeader.iloc[0,3]
tFin=dfHeader.iloc[0,4]
numIntervalosT=dfHeader.iloc[0,5]
#M00=dfHeader.iloc[0,6]

#g=pandas.read_csv(datosTxt, sep=';', skiprows=2, header=None, chunksize=100, iterator=True)

g=pandas.read_csv(datosTxt, sep=';', skiprows=2, error_bad_lines=False)
g.info(memory_usage="deep")
g=g.to_numpy()
g=g.astype(np.float) 

dfX12=pandas.read_csv(datosTxt, skiprows=1, nrows=2, header=None, sep=';')

xi12=dfX12.iloc[0,:]
xi12=xi12.to_numpy()

'''
#Recurso de uso caso de que se quiera comparar con una solucion exacta
xi_exacta=np.zeros(20000)
xInicial=xi12[0]
xFinal=20
for i in range(0, len(xi_exacta)):
    xi_exacta[i]=xInicial + i * (xFinal-xInicial)/(len(xi_exacta)-1)
'''

xi=np.zeros(len(xi12)-1)
for i in range(0,len(xi)):
    xi[i]=0.5*(xi12[i]+xi12[i+1])


##############
file1 = open(datosTxt)
lastLines1 = tl.tail(file1,3)
file1.close()
df1=pandas.read_csv(io.StringIO('\n'.join(lastLines1)), sep=';', header=None, error_bad_lines=False)

tiempo=df1.iloc[1]
tiempo=tiempo.to_numpy()

#En segundos
tiempo1=0.25
tiempo2=0.8
tiempo3=1
tiempo4=1.75


posiciones=[]
#print(tiempo)
#Busca el instante mas cercano en el array de tiempos generado en el fichero de texto de la simulacion
for i in range(0,len(tiempo)):
    if(tiempo1>=tiempo[i] and tiempo1<tiempo[i+1]):
        pos1=i
        posiciones.append(pos1)
        print('Tiempo 1: ', tiempo[i])
    if(tiempo2>=tiempo[i] and tiempo2<tiempo[i+1]):
        pos2 = i
        posiciones.append(pos2)
        print('Tiempo 2: ', tiempo[i])
    if(tiempo3>=tiempo[i] and tiempo3<tiempo[i+1]):
        pos3 = i
        posiciones.append(pos3)
        print('Tiempo 3: ', tiempo[i])
    if (tiempo4 >= tiempo[i] and tiempo4 < tiempo[i+1]):
        pos4 = i
        posiciones.append(pos4)
        print('Tiempo 4: ', tiempo[i])



plt.figure(1)

plt.plot(xi,xi*g[pos1,:],'b',label="t="+str(tiempo1))
plt.plot(xi,xi*g[pos2,:],'g',label="t="+str(tiempo2))
plt.plot(xi,xi*g[pos3,:],'r',label="t="+str(tiempo3))
plt.plot(xi,xi*g[pos4,:],'m',label="t="+str(tiempo4))
plt.xlabel('x')
plt.ylabel('xn(x,t)')
plt.legend()
plt.title('Seccional (Ker(x,y)=1)')
#print(pos1, pos2, pos3, pos4)

plt.figure(2)
plt.plot(np.log10(xi),np.log10(xi*g[pos1,:]),'b',label="t="+str(tiempo1))
####
plt.plot(np.log10(xi),np.log10(xi*g[pos2,:]),'g',label="t="+str(tiempo2))
####
plt.plot(np.log10(xi),np.log10(xi*g[pos3,:]),'r',label="t="+str(tiempo3))
####
plt.plot(np.log10(xi),np.log10(xi*g[pos4,:]),'m',label="t="+str(tiempo4))


plt.legend()
plt.title('Seccional (Ker(x,y)=1)')
plt.xlabel('log(x)')
plt.ylabel('log(xn(x,t))')
#plt.ylabel('n(x,t)')
#plt.xlim(-0.2,12.2)
#plt.xlim(0,250)
#plt.ylim(-0.001,1.2)
#plt.ylim(0,0.8)

plt.show()

