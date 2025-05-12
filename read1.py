import os
import matplotlib.pyplot as plt


eV_data = []
Jclas_data = []
Prob3clas = []
Prob4clas = []
Prob5clas = []
Prob6clas = []
Prob7clas = []
Prob8clas = []
Iclas_data = []


Jqm_data = []
Prob3qm = []
Prob4qm = []
Prob5qm = []
Prob6qm = []
Prob7qm = []
Prob8qm = []
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
    #print(float(item))  
    aux0.append(0)
nuevo.close()
        
nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    Prob3clas.append(float(item))  
    #print(float(item)) 
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
    Prob7clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[7] for data in nuevo]:  
    Prob8clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[8] for data in nuevo]:  
    Iclas_data.append(float(item))   
nuevo.close()


datos0 = 'lindbladgamU'
nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    Jqm_data.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)     
for item in [data.split()[2] for data in nuevo0]:  
    Prob3qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[3] for data in nuevo0]:  
    Prob4qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[4] for data in nuevo0]:  
    Prob5qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[5] for data in nuevo0]:  
    Prob6qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[6] for data in nuevo0]:  
    Prob7qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[7] for data in nuevo0]:  
    Prob8qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[8] for data in nuevo0]:  
    Iqm_data.append(float(item))   
nuevo0.close()

plt.plot(eV_data,Jclas_data, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Jqm_data, color='blue',lw = 4, label = "quantum")
plt.plot(eV_data,aux0)
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()


plt.plot(eV_data,Prob3clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob3qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{101}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob4clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob4qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{100}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob5clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob5qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{011}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob6clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob6qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{010}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob7clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob7qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{001}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob8clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob8qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{000}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Iclas_data, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Iqm_data, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{D}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
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