import os
import matplotlib.pyplot as plt

eV_data = []
Nparcial = []
Jlparcial = []
Jrparcial = []
JlRparcial = []
coheparcial = []
concuparcial = []

Nred = []
Jlred = []
Jrred = []
JlRred = []
cohered = []
concured = []

#######################################
########compararcodigos################
#datos = 'parcialweakg'
datos = 'parcialstrongg'
#datos = 'parcialweakg'
fichero = open(datos)
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    eV_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    Nparcial.append(-float(item)) 
    #print(float(item))  
nuevo.close()

nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    Jlparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()

nuevo = open(datos)
for item in [data.split()[3] for data in nuevo]:  
    Jrparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()
        
nuevo = open(datos)
for item in [data.split()[4] for data in nuevo]:  
    JlRparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()
nuevo = open(datos)
for item in [data.split()[5] for data in nuevo]:  
    coheparcial.append(float(item))   
nuevo.close()
nuevo = open(datos)
for item in [data.split()[6] for data in nuevo]:  
    concuparcial.append(float(item))   
nuevo.close()

#datos0 = 'redfieldweakg'
datos0 = 'redfieldstrongg'
#datos0 = 'redfieldweakg'
fichero0 = open(datos0)

for item in [data.split()[1] for data in fichero0]: #saca la columna 0|
    Nred.append(-float(item))    
fichero0.close()
fichero0 = open(datos0)
for item in [data.split()[2] for data in fichero0]: #saca la columna 0|
    Jlred.append(float(item))
fichero0.close()    

fichero0 = open(datos0) 
for item in [data.split()[3] for data in fichero0]: #saca la columna 0|
    Jrred.append(float(item))
fichero0.close()

fichero0 = open(datos0)
for item in [data.split()[4] for data in fichero0]: #saca la columna 0|
    JlRred.append(float(item))      
fichero0.close()
fichero0 = open(datos0)
for item in [data.split()[5] for data in fichero0]: #saca la columna 0|
    cohered.append(float(item))
fichero0.close()
fichero0 = open(datos0)
for item in [data.split()[6] for data in fichero0]: #saca la columna 0|
    concured.append(float(item))
fichero0.close()        





plt.plot(eV_data,Nparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Nred, color='blue',lw = 4, label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jlparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Jlred, color='blue',lw = 4, label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jrparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Jrred, color='blue',lw = 4, label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{R}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,JlRparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,JlRred, color='blue',lw = 4, label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{LR}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,coheparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,cohered, color='blue',lw = 4,  label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$C$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,concuparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,concured, color='blue',lw = 4, label = "redfield")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$Concu$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()
