import numpy as np
import matplotlib.pyplot as plt
import pandas
from matplotlib import colors
import tailer as tl

import io

#Documento txt a representar:
#datosTxt='volumentesFin.txt'
#datos1Txt='26JUN/ID/SECC_ID_m150Equi_min5int45RK3.txt'

datos1Txt='FEBRERO/vol_finitos/trapecios_comp/RK3/FINAL/2s-vol_finitos_RK3_SUMA_trapecios_compuesto_50.txt'
datos2Txt='FEBRERO/vol_finitos/trapecios_comp/RK3/FINAL/2s-vol_finitos_RK3_SUMA_trapecios_compuesto_100.txt'
datos3Txt='FEBRERO/vol_finitos/trapecios_comp/RK3/FINAL/2s-vol_finitos_RK3_SUMA_trapecios_compuesto_200.txt'




dfHeader = pandas.read_csv(datos1Txt, sep=';',delimiter=';',nrows=1,names=["0","1","2","3","4","5","6"])

t0=dfHeader.iloc[0,3]
tFin=dfHeader.iloc[0,4]

##############
file1 = open(datos1Txt)
lastLines1 = tl.tail(file1,3) #to read last 15 lines, change it  to any value.
file1.close()
df1=pandas.read_csv(io.StringIO('\n'.join(lastLines1)), sep=';', header=None, error_bad_lines=False)

t_Q1=df1.iloc[1]
tQ1=t_Q1.to_numpy()
#print(tQ1[11000])


Q1 = pandas.read_csv(datos1Txt, skiprows=4, sep=';', header=None, error_bad_lines=False)

TQ_1=np.zeros(len(tQ1))
errores_1 = df1.iloc[0]
errores_1=errores_1.to_numpy()
#print(errores_1)
errores1=np.zeros(len(tQ1))
for i in range(0, len(tQ1)):
    errores1[i]=errores_1[i]
    TQ_1[i]=tQ1[i]


##############
file2 = open(datos2Txt)
lastLines2 = tl.tail(file2,3) #to read last 15 lines, change it  to any value.
file2.close()
df2=pandas.read_csv(io.StringIO('\n'.join(lastLines2)), sep=';', header=None, error_bad_lines=False)

t_Q2=df2.iloc[1]
tQ2=t_Q2.to_numpy()

Q2 = pandas.read_csv(datos2Txt, skiprows=4, sep=';', header=None, error_bad_lines=False)

TQ_2=np.zeros(len(tQ2))
errores_2 = df2.iloc[0]
errores_2=errores_2.to_numpy()
errores2=np.zeros(len(tQ2))
for i in range(0, len(tQ2)):
    errores2[i]=errores_2[i]
    TQ_2[i] = tQ2[i]



#################

file3 = open(datos3Txt)
lastLines3 = tl.tail(file3,3) #to read last 15 lines, change it  to any value.
file3.close()
df3=pandas.read_csv(io.StringIO('\n'.join(lastLines3)), sep=';', header=None, error_bad_lines=False)

t_Q3=df3.iloc[1]
tQ3=t_Q3.to_numpy()

Q3 = pandas.read_csv(datos3Txt, skiprows=4, sep=';', header=None, error_bad_lines=False)

TQ_3=np.zeros(len(tQ3))
errores_3 = df3.iloc[0]
errores_3=errores_3.to_numpy()
errores3=np.zeros(len(tQ3))
for i in range(0, len(tQ3)):
    errores3[i]=errores_3[i]
    TQ_3[i] = tQ3[i]
    if(errores3[i]>11):
        print("Posicion del error", i)


#################


plt.figure(1)

plt.plot(tQ1,np.log10(errores_1), 'b', label='Nx=50')
plt.plot(tQ2,np.log10(errores_2), 'g', label='Nx=100')
plt.plot(tQ3,np.log10(errores_3), 'r', label='Nx=200')
plt.xlabel('t')
plt.title('Errores volúmenes finitos (ker(x,y)=x+y)')
plt.legend()


plt.figure(2)

plt.plot(tQ1,errores_1, 'b', label='Nx=50')
plt.plot(tQ2,errores_2, 'g', label='Nx=100')
plt.plot(tQ3,errores_3,'r', label='Nx=200')
#plt.xlim(0,0.2)
#plt.ylim(0,1.5)
plt.xlabel('t')
plt.title('Errores volúmenes finitos (ker(x,y)=x+y)')
#plt.xlim(0.8,1.5)
#plt.ylim(0.3,0.6)
plt.legend()
plt.show()

def meanError(x):
    return np.mean(x)

def stddevError(x):
    return np.std(x)

def maxError(x):
    return np.max(x)

print("1------------")
print("Max=",maxError(errores1))
print("Mean=",meanError(errores1))
print("Stddev=",stddevError(errores1))

print("2------------")
print("Max=",maxError(errores2))
print("Mean=",meanError(errores2))
print("Stddev=",stddevError(errores2))

print("3------------")
print("Max=",maxError(errores3))
print("Mean=",meanError(errores3))
print("Stddev=",stddevError(errores3))