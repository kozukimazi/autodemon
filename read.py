import os
import matplotlib.pyplot as plt


eV_data = []
Jclas_data = []
Iclas_data = []

Jqm_data = []
Iqm_data = []

datos = 'stochastic'
fichero = open("stochastic")
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    eV_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    Jclas_data.append(-float(item))   
nuevo.close()
        
nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    Iclas_data.append(float(item))   
nuevo.close()

datos0 = 'lindblad'
nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    Jqm_data.append(float(item))   
nuevo0.close()
        
nuevo0 = open(datos0)
for item in [data.split()[2] for data in nuevo0]:  
    Iqm_data.append(float(item))   
nuevo0.close()

plt.scatter(eV_data,Jclas_data,s = 1, label = r'$\mathcal{J}_{cl}$')
plt.scatter(eV_data,Jqm_data,s = 1, label = r'$\mathcal{J}_{qm}$')
plt.legend()
plt.show()

plt.scatter(eV_data,Iclas_data,s = 0.5, label = r'$\mathcal{I}_{cl}$')
plt.scatter(eV_data,Iqm_data,s=0.5, label = r'$\mathcal{I}_{qm}$')
plt.legend()
plt.show()

Num = len(eV_data)
difI = []
difJ = []
for i in range(Num):
    aux0 =Jclas_data[i] - Jqm_data[i] 
    aux1 = Iclas_data[i] - Iqm_data[i]
    difJ.append(aux0)
    difI.append(aux1)

plt.plot(eV_data,difI, label = r'$\mathcal{I}_{cl} - \mathcal{I}_{qm}$')
plt.legend()
plt.show()

plt.plot(eV_data,difJ, label = r'$\mathcal{J}_{cl} - \mathcal{J}_{qm}$')
plt.legend()
plt.show()
