import numpy as np
import matplotlib.pyplot as plt
import pandas
import scipy.integrate as integrar
import tailer as tl
import io
from decimal import Decimal

###################################
### Representacion grafica de las aproximaciones obtenidas en los casos de Oceanografia
### Figuras del capitulo 3 de la memoria del trabajo
###################################


datosTxt='ficheros_simulaciones/OCEANOGRAFIA/vol_finitos/VOL-200_e-13mol_Z3_6m_EXPONENCIAL.txt'


dfHeader = pandas.read_csv(datosTxt, sep=';',delimiter=';',nrows=1,names=["0","1","2","3","4","5","6"])

Nx=dfHeader.iloc[0,0]
iniX12=dfHeader.iloc[0,1]
R=dfHeader.iloc[0,2]
t0=dfHeader.iloc[0,3]
tFin=dfHeader.iloc[0,4]
numIntervalosT=dfHeader.iloc[0,5]
#M00=dfHeader.iloc[0,6]

columnas=Nx+1

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
tiempo1=200000
tiempo2=1000000
tiempo3=8000000
tiempo4=12000000


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


add0=0

#Funcion para convertir los resultados de las simulaciones de la escala en masa a diametro para su posterior graficado
def MasaToDiam(x, g, posicion):
    D=3.0
    m0=x[0]
    p0=1.17e11
    d0=2*(3*m0/(4*np.pi*p0))**(1./3)
    xDiam=np.zeros(len(x))
    nDiam=np.zeros(len(x))
    for i in range(0,len(x)):
        xDiam[i]=d0*(x[i]/m0)**(1./D)
        if(Decimal(g[posicion,i])<1e-60):
            g[posicion, i]=0
        try:
            nDiam[i]=Decimal(m0*D/d0*(xDiam[i]/d0)**(D-1))*Decimal(g[posicion,i]+add0)
            nDiam[i]=(nDiam[i]/(6.02214129*10**23))*(10**6)*(10**(-3))
        except:
            print(i, ': ',g[posicion,i], ', ', Decimal(m0*D/d0*(xDiam[i]/d0)**(D-1))*Decimal(g[posicion,i]))
    return xDiam,nDiam


plt.figure(1)


tiempo1=int(tiempo1/86400)
tiempo2=int(tiempo2/86400)
tiempo3=int(tiempo3/86400)
tiempo4=int(tiempo4/86400)


plt.plot(xi,xi*g[pos1,:],'b',label="t="+str(tiempo1)+" d")
plt.plot(xi,xi*g[pos2,:],'g',label="t="+str(tiempo2)+" d")
plt.plot(xi,xi*g[pos3,:],'r',label="t="+str(tiempo3)+" d")
plt.plot(xi,xi*g[pos4,:],'m',label="t="+str(tiempo4)+" d")

print(pos1, pos2, pos3, pos4)

plt.figure(2)
####
plot1=MasaToDiam(xi,g,pos1)
plt.plot(10**6*plot1[0],plot1[1],'b',label="t="+str(tiempo1)+" d")
####
plot2=MasaToDiam(xi,g,pos2)
plt.plot(10**6*plot2[0],plot2[1],'g',label="t="+str(tiempo2)+" d")
####
plot3=MasaToDiam(xi,g,pos3)
plt.plot(10**6*plot3[0],plot3[1],'r',label="t="+str(tiempo3)+" d")
####
plot4=MasaToDiam(xi,g,pos4)
plt.plot(10**6*plot4[0],plot4[1],'m',label="t="+str(tiempo4)+" d")
#####
plt.figure(3)
add1=0

plt.plot(np.log10(10**6*plot1[0]),np.log10(np.abs(plot1[1]+add1)),'b',label="t="+str(tiempo1)+" d")
####
plt.plot(np.log10(10**6*plot2[0]),np.log10(np.abs(plot2[1]+add1)),'g',label="t="+str(tiempo2)+" d")
####
plt.plot(np.log10(10**6*plot3[0]),np.log10(plot3[1]+add1),'r',label="t="+str(tiempo3)+" d")
####
plt.plot(np.log10(10**6*plot4[0]),np.log10(plot4[1]+add1),'m',label="t="+str(tiempo4)+" d")

'''
plt.plot(10**6*plot1[0],np.log10(plot1[1]+add1),'b',label="t="+str(tiempo1)+" d")
####
plt.plot(10**6*plot2[0],np.log10(plot2[1]+add1),'g',label="t="+str(tiempo2)+" d")
####
plt.plot(10**6*plot3[0],np.log10(plot3[1]+add1),'r',label="t="+str(tiempo3)+" d")
####
plt.plot(10**6*plot4[0],np.log10(plot4[1]+add1),'m',label="t="+str(tiempo4)+" d")
'''


#plt.plot(xi,g[pos2,:],'g',label="t="+str(tiempo2))
#plt.plot(xi,g[pos3,:],'r',label="t="+str(tiempo3))
plt.legend()
plt.title('Z=3 m, $\mathcal{I}$=0.1 $mg L^{-1} d^{-1}$')
plt.xlabel('d (diámetro, en $\mu m)$')
plt.ylabel('$\mu$moles $L^{-1}$')
#plt.ylabel('n(x,t)')
#plt.xlim(-0.2,12.2)
#plt.xlim(0,250)
#plt.ylim(-0.001,1.2)
#plt.ylim(0,0.8)

#for i in range(0,len(g[2*pos1,:])):
    #if(g[2*pos1+1, i]==0):
     #   print(g[2*pos1,i])


plt.figure(4)
add2=0

plt.plot(np.log10(xi),np.log10(g[pos1,:]+add2),'b',label="t="+str(tiempo1)+" d")
plt.plot(np.log10(xi),np.log10(g[pos2,:]+add2),'g',label="t="+str(tiempo2)+" d")
plt.plot(np.log10(xi),np.log10(g[pos3,:]+add2),'r',label="t="+str(tiempo3)+" d")
plt.plot(np.log10(xi),np.log10(g[pos4,:]+add2),'m',label="t="+str(tiempo4)+" d")

plt.legend()
plt.title('Z=30 m, $\mathcal{I}$=0.01 $mg l^{-1} d^{-1}$')
plt.xlabel('log(x) (masa, en $\mu g)$')
plt.ylabel('log(n(x,t))')

'''
plt.figure(2)
#plt.plot(np.log(xi_exacta),np.log(yG_ini),'k', label='Exacta (t=0)')
plt.plot(np.log(xi),np.log(g[pos1,:]),'b',label="t="+str(tiempo1))
plt.plot(np.log(xi),np.log(g[pos2,:]),'g',label="t="+str(tiempo2))
plt.plot(np.log(xi),np.log(g[pos3,:]),'r',label="t="+str(tiempo3))
plt.legend()
plt.title('Volúmenes finitos (ker(x,y)=$(x^{1/3}+y^{1/3})^2|x^{1/3}-y^{1/3}|$, Nx=300)')
plt.xlabel('log(x)')
plt.ylabel('log(x u(x))')
#plt.xlim(-5,-3)
#plt.ylim(-5,1)
'''

plt.show()

