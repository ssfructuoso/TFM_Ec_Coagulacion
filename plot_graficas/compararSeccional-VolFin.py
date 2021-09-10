import pandas
import numpy
import matplotlib.pyplot as plt

###########################################
#### Script utilizado para comparar la diferencia de los resultados obtenidos con el metodo de seccional frente al método
#### de volumenes finitos en las simulaciones mostradas al final del capítulo 3 para el caso de Oceanografia y con NX=50.
#### Los resultados se muestran en la figura 3.5 de la memoria del trabajo.
###########################################


datosSecc='ficheros_simulaciones/OCEANOGRAFIA/comparativa_metodos/SECC_test-100_e-13mol_Z3_6m_POTENCIAS_1m.txt'
datosVolFin='ficheros_simulaciones/OCEANOGRAFIA/comparativa_metodos/VOL_TESTs_e-13mol_Z3_6m_POTENCIAS_1m.txt'


data_Secc = pandas.read_csv(datosSecc, sep=';',delimiter=';', skiprows=2, skipfooter=3, header=None)
data_Vol = pandas.read_csv(datosVolFin, sep=';',delimiter=';', skiprows=2, skipfooter=3, header=None)

data_Secc.to_numpy()
data_Vol.to_numpy()

#1 mes
time=[0,300,600,900,1200,1500,1800,2100,2400,2700,3000,3300,3600,3900,4200,4500,4800,5100,5400,5700,6000,6300,6600,6900,7200,7500,7800,8100,8400,8700,9000,9300,9600,9900,10200,10500,10800,11100,11400,11700,12000,12300,12600,12900,13200,13500,13800,14100,14400,14700,15000,15300,15600,15900,16200,16500,16800,17100,17400,17700,18000,18300,18600,18900,19200,19500,19800,20100,20400,20700,21000,21300,21600,21900,22200,22500,22800,23100,23400,23700,24000,24300,24600,24900,25200,25500,25800,26100,26400,26700,27000,27300,27600,27900,28200,28500,42900,57300,71700,86100,100500,114900,129300,143700,158100,172500,186900,201300,215700,230100,244500,258900,273300,287700,302100,316500,330900,345300,359700,374100,388500,402900,417300,431700,446100,460500,474900,489300,503700,518100,532500,546900,561300,575700,590100,604500,618900,633300,647700,662100,676500,690900,705300,719700,734100,748500,762900,777300,791700,806100,820500,834900,849300,863700,878100,892500,906900,921300,935700,950100,964500,978900,993300,1007700,1022100,1036500,1050900,1065300,1079700,1094100,1108500,1122900,1137300,1151700,1166100,1180500,1194900,1209300,1223700,1238100,1252500,1266900,1281300,1295700,1310100,1324500,1338900,1353300,1367700,1382100,1396500,1410900,1425300,1439700,1454100,1468500,1482900,1497300,1511700,1526100,1540500,1554900,1569300,1583700,1598100,1612500,1626900,1641300,1655700,1670100,1684500,1698900,1713300,1727700,1742100,1756500,1770900,1785300,1799700,1814100,1828500,1842900,1857300,1871700,1886100,1900500,1914900,1929300,1943700,1958100,1972500,1986900,2001300,2015700,2030100,2044500,2058900,2073300,2087700,2102100,2116500,2130900,2145300,2159700,2174100,2188500,2202900,2217300,2231700,2246100,2260500,2274900,2289300,2303700,2318100,2332500,2346900,2361300,2375700,2390100,2404500,2418900,2433300,2447700,2462100,2476500,2490900,2505300,2519700,2534100,2548500,2562900,2577300]

max=numpy.zeros(len(time))

#print(data_Vol)
#print(data_Secc)
n=50
diff=numpy.zeros(n)
for i in range(0,len(max)):
    for j in range(0,n):
        try:
            diff[j]=numpy.abs(100*(data_Secc[j][i]-data_Vol[j][i])/((data_Vol[j][i])))
        except:
            print('Error: ',data_Secc[j][i],',',data_Vol[j][i])

    max[i]=numpy.max(diff)

#print(max)


for i in range(0,len(time)):
    time[i]=time[i]/86400

plt.figure(1)

plt.plot(time,max)
plt.title('n(x,t$_{0}$)=$6.023\cdot 10^{10}$(x$_{0}$/x)$^{1/5}$ POC/m$^{3}$')
plt.ylabel('R(t)')
plt.xlabel('t (días)')

plt.show()