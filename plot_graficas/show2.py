import numpy as np
import matplotlib.pyplot as plt
import pandas
import math
import scipy.integrate as integrar
import tailer as tl
import io
from decimal import Decimal
import csv
from matplotlib import colors


#Documento txt a representar:
#datosTxt='volumentesFin.txt'
datosTxt='FEBRERO/seccional/RK3/midpoint/seccional_RK3_SUMA-mid_point_200.txt'


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
#getData = csv.reader(open(datosTxt))

#chunk, chunksize = [], 100


dfX12=pandas.read_csv(datosTxt, skiprows=1, nrows=2, header=None, sep=';')


xi12=dfX12.iloc[0,:]
xi12=xi12.to_numpy()
#print(xi12)

xi_exacta=np.zeros(20000)
xInicial=xi12[0]
xFinal=20
for i in range(0, len(xi_exacta)):
    xi_exacta[i]=xInicial + i * (xFinal-xInicial)/(len(xi_exacta)-1)

xi=np.zeros(len(xi12)-1)
for i in range(0,len(xi)):
    xi[i]=0.5*(xi12[i]+xi12[i+1])



#xi[len(xi)-1]=xi12[len(xi12)-1]#Revisar!!! imprecisoooo
#my_cols = [str(i) for i in range(filas*columnas)]

#df = pandas.read_csv(datosTxt, sep=";", names=my_cols)
#g=g.to_numpy();
#g=np.delete(g,0,axis=0)


##############
file1 = open(datosTxt)
lastLines1 = tl.tail(file1,3)
file1.close()
df1=pandas.read_csv(io.StringIO('\n'.join(lastLines1)), sep=';', header=None, error_bad_lines=False)

tiempo=df1.iloc[1]
tiempo=tiempo.to_numpy()

#print(tiempo)

#############
def I_bessel(x):
    integrando=lambda theta: np.exp(x*np.cos(theta))*np.cos(theta)
    integral=integrar.quad(integrando, 0, np.pi)[0]
    integral=integral/(np.pi)
    return integral

#print(I_bessel(0))
'''
def gIni(x,t,M00):
    #return x*((2*M00)/(2+M00*t))**2 *np.exp(-2*M00*x/(2+M00*t))
    #return x*(2*math.pi)**(-1.0/2)*np.exp(-t)*x**(-3.0/2)*np.exp(-x/2*np.exp(-2*t))

    T=0
    if(t<=1):
        T=1+t
    else:
        T=2*np.sqrt(t)

    return x * np.exp(-T*x)*I_bessel(2*x*np.sqrt(t))/(x**2 * np.sqrt(t))
####################
'''

tiempo1=0
tiempo2=0.187
tiempo3=0.5
tiempo4=1

posiciones=[]
pos1=0
pos2=0
pos3=0
pos4=0
print(tiempo)
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

'''
yData=[]
for i, line in enumerate(getData):
    if (i % chunksize == 0 and i > 0):
        for p in posiciones:
           if(p>=i and p<i+chunksize):
               print(chunk)
        del chunk[:]  # or: chunk = []
    chunk.append(line)
'''

def MasaToDiam(x, g, posicion):
    D=2.6
    m0=x[0]
    p0=2.25e12
    d0=2*(3*m0/(4*np.pi*p0))**(1./3)
    xDiam=np.zeros(len(x))
    nDiam=np.zeros(len(x))
    for i in range(0,len(x)):
        xDiam[i]=d0*(x[i]/m0)**(1./D)
        if(Decimal(g[posicion,i])<1e-60):
            g[posicion, i]=0
        try:
            nDiam[i]=Decimal(m0*D/d0*(xDiam[i]/d0)**(D-1))*Decimal(g[posicion,i])
            nDiam[i]=(nDiam[i]/(6.02214129*10**23))*(10**6)*(10**(-3))
        except:
            print(i, ': ',g[posicion,i], ', ', Decimal(m0*D/d0*(xDiam[i]/d0)**(D-1))*Decimal(g[posicion,i]))
    return xDiam,nDiam


plt.figure(1)
'''
yG_ini=np.zeros(len(xi_exacta))
yG_ini2=np.zeros(len(xi_exacta))
for i in range(0,len(xi_exacta)):
    yG_ini[i]=gIni(xi_exacta[i],tiempo[pos1],M00)
    #yG_ini2[i] = gIni(xi_exacta[i], tiempo[pos3], M00)
'''

#plt.plot(xi_exacta,yG_ini,'k', label='Exacta (t=0)')
#plt.plot(xi_exacta,yG_ini2,'g', label='Exacta (t=0)')
#print(g[2*pos1,:])
#print(g[2*pos1+1,:])
'''
tiempo1=int(tiempo1/60)
tiempo2=int(tiempo2/60)
tiempo3=int(tiempo3/60)
tiempo4=int(tiempo4/60)
'''
pos1=2
pos2=10
pos3=20
pos4=150

plt.plot(xi,g[pos1,:],'b',label="t="+str(tiempo1)+" min")
plt.plot(xi,g[pos2,:],'g',label="t="+str(tiempo2)+" min")
plt.plot(xi,g[pos3,:],'r',label="t="+str(tiempo3)+" min")
plt.plot(xi,g[pos4,:],'m',label="t="+str(tiempo4)+" min")

print(pos1, pos2, pos3, pos4)

plt.figure(2)
####
plot1=MasaToDiam(xi,g,pos1)
plt.plot(10**6*plot1[0],plot1[1],'b',label="t="+str(tiempo1)+" min")
####
plot2=MasaToDiam(xi,g,pos2)
plt.plot(10**6*plot2[0],plot2[1],'g',label="t="+str(tiempo2)+" min")
####
plot3=MasaToDiam(xi,g,pos3)
plt.plot(10**6*plot3[0],plot3[1],'r',label="t="+str(tiempo3)+" min")
####
plot4=MasaToDiam(xi,g,pos4)
plt.plot(10**6*plot4[0],plot4[1],'m',label="t="+str(tiempo4)+" min")
#####


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


plt.figure(3)

plt.plot(np.log(xi),np.log(g[pos1,:]),'b',label="t="+str(tiempo1)+" min")
plt.plot(np.log(xi),np.log(g[pos2,:]),'g',label="t="+str(tiempo2)+" min")
plt.plot(np.log(xi),np.log(g[pos3,:]),'r',label="t="+str(tiempo3)+" min")
plt.plot(np.log(xi),np.log(g[pos4,:]),'m',label="t="+str(tiempo4)+" min")

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

