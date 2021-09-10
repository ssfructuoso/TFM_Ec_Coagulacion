import numpy as np
import matplotlib.pyplot as plt
import pandas
import scipy.integrate as integrar
import tailer as tl
import io
from scipy.interpolate import CubicSpline


###########################################
### Script utilizado para calcular el primer momento M1(t) con los dos metodos numericos desarrollados.
### En este caso se ha escogido el nucleo producto, pudiendose diferenciar el fenomeno de gelacion.
### Este script es utilizado para generar la figura 2.10 de la memoria del trabajo.
###########################################


def plot1(datosTxt, pltlabel):
    
    dfHeader = pandas.read_csv(datosTxt, sep=';', delimiter=';', nrows=1, names=["0", "1", "2", "3", "4", "5", "6"])

    iniX12 = dfHeader.iloc[0, 1]
    R = dfHeader.iloc[0, 2]

    g = pandas.read_csv(datosTxt, skiprows=2, sep=';', header=None, error_bad_lines=False)

    dfX12 = pandas.read_csv(datosTxt, skiprows=1, nrows=2, header=None, sep=';', error_bad_lines=False)

    xi12 = dfX12.iloc[0, :]
    xi12 = xi12.to_numpy()

    xi = np.zeros(len(xi12) - 1)
    for i in range(0, len(xi)):
        xi[i]=0.5*(xi12[i]+xi12[i+1])

    g = g.to_numpy()

    ##############
    file1 = open(datosTxt)
    lastLines1 = tl.tail(file1, 3)
    file1.close()
    df1 = pandas.read_csv(io.StringIO('\n'.join(lastLines1)), sep=';', header=None, error_bad_lines=False)

    tiempo = df1.iloc[1]
    tiempo = tiempo.to_numpy()

    #############
    def m1(posicion):
        x=np.append(xi12[0],xi)
        x=np.append(x,xi12[len(xi12)-1])
        y = np.append(g[posicion, 0], g[posicion, :])
        y = np.append(y, g[posicion , len(g[posicion, :]) - 1])
        # f = interp1d(x, y, kind='cubic')
        f = CubicSpline(x, x*y)
        # f=CubicSpline(x,y)
        # f_integrando = lambda xx: f(xx)
        # integral_f=integrar.quad(f_integrando, iniX12, iniX12+ R)[0]
        x_simp = np.linspace(iniX12, iniX12 + R, 100000)
        integral_f = integrar.simps(f(x_simp), x_simp)
        return integral_f

    fTiempo = np.zeros(len(tiempo))
    for i in range(0, len(tiempo)):
        fTiempo[i] = m1(i)

    plt.figure(1)
    plt.plot(tiempo, fTiempo, label=pltlabel)
    
######################
######################

#### Momento M1 exacto para el nucleo producto (Ker(x,y)=xy)
tExacta=np.linspace(0.001,2,num=10000)
m1Exacta=[]
for i in tExacta:
    if i<=1:
        m1Exacta.append(1)
    else:
        m1Exacta.append(i**(-0.5))


plt.plot(tExacta,m1Exacta, 'k', label='Sol. exacta')###
plt.title('M1 (Ker(x,y)=xy)')
plt.xlabel('t')
plt.ylabel('M$_{1}$(t)')
plt.xlim(0,1.6)
plt.ylim(0.6,1.2)


datosVol= 'ficheros_simulaciones/M1/Ker_producto/vol_finitos_RK3_PROD_trapecios_compuesto_100.txt'
datosSec='ficheros_simulaciones/M1/Ker_producto/secc_RK3_PROD_trapecios_compuesto_100.txt'
plot1(datosVol, "Vol. finitos") #Momento M1 para el nucleo producto con el metodo de volumenes finitos
plot1(datosSec, "Seccional")#Momento M1 para el nucleo producto con el metodo seccional
xV=[1,1]
yV=[-1,2]
plt.plot(xV,yV,'r--')
plt.legend()
#plt.gca().xaxis.set_major_formatter(Formatter(fformat="%1.1f"))
plt.xlim(0,2.5)
plt.show()
