import os
import matplotlib.pyplot as plt


eV_data = []
Jclas_data = []
Prob3clas = []
Prob4clas = []
Prob5clas = []
Prob6clas = []
Iclas_data = []


Jqm_data = []
Prob3qm = []
Prob4qm = []
Prob5qm = []
Prob6qm = []
Iqm_data = []
aux0 = []
#aqui se abre el archivo de datos que quiero usar para comparar
datos = 'classic'
fichero = open("classic")
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    eV_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    Jclas_data.append(-float(item))   
    aux0.append(0)
nuevo.close()
        
nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    Prob3clas.append(float(item))   
nuevo.close()
nuevo = open(datos)
for item in [data.split()[3] for data in nuevo]:  
    Prob4clas.append(float(item))   
nuevo.close()
nuevo = open(datos)
for item in [data.split()[4] for data in nuevo]:  
    Prob5clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[5] for data in nuevo]:  
    Prob6clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[6] for data in nuevo]:  
    Iclas_data.append(float(item))   
nuevo.close()

datos0 = 'lindblad'
nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    Jqm_data.append(float(item))   
nuevo0.close()
        
for item in [data.split()[2] for data in nuevo0]:  
    Prob3qm.append(float(item))   
nuevo0.close()

for item in [data.split()[3] for data in nuevo0]:  
    Prob4qm.append(float(item))   
nuevo0.close()

for item in [data.split()[4] for data in nuevo0]:  
    Prob5qm.append(float(item))   
nuevo0.close()

for item in [data.split()[5] for data in nuevo0]:  
    Prob6qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[6] for data in nuevo0]:  
    Iqm_data.append(float(item))   
nuevo0.close()

plt.plot(eV_data,Jclas_data, color='red',lw = 4, label = r'$\mathcal{J}_{cl}$')
plt.plot(eV_data,Jqm_data, color='blue',lw = 4, label = r'$\mathcal{J}_{qm}$')
plt.plot(eV_data,aux0)
plt.legend()
plt.show()

plt.plot(eV_data,Prob3clas, color='red',lw = 4, label = r'$\rho_{101}_{cl}$')
plt.plot(eV_data,Prob3qm, color='blue',lw = 4, label = r'$\\rho_{101}_{qm}$')
plt.legend()
plt.show()

plt.plot(eV_data,Prob4clas, color='red',lw = 4, label = r'$\rho_{100}_{cl}$')
plt.plot(eV_data,Prob4qm, color='blue',lw = 4, label = r'$\rho_{100}_{qm}$')
plt.legend()
plt.show()

plt.plot(eV_data,Prob5clas, color='red',lw = 4, label = r'$\rho_{011}_{cl}$')
plt.plot(eV_data,Prob5qm, color='blue',lw = 4, label = r'$\rho_{011}_{qm}$')
plt.legend()
plt.show()

plt.plot(eV_data,Prob6clas, color='red',lw = 4, label = r'$\rho_{010}_{cl}$')
plt.plot(eV_data,Prob6qm, color='blue',lw = 4, label = r'$\rho_{010}_{qm}$')
plt.legend()
plt.show()

plt.plot(eV_data,Iclas_data, color='red',lw = 4, label = r'$\mathcal{I}_{cl}$')
plt.plot(eV_data,Iqm_data, color='blue',lw = 4, label = r'$\mathcal{I}_{qm}$')
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
