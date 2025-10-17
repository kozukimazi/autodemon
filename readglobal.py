import os
import matplotlib.pyplot as plt

g_data = []
coheparcial = []
concuparcial = []


cohered = []
concured = []


coheglobal = []
concuglobal = []

#######################################
########compararcodigos################
#datos = 'parcialweakg'
datos = 'semicouplarge'
#datos = 'parcialweakg'
fichero = open(datos)
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    g_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    concuparcial.append(float(item)) 
    #print(float(item))  
nuevo.close()

nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    coheparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()


#######################################
########compararcodigos################
#datos = 'parcialweakg'
datos0 = 'redcouplarge'
#datos = 'parcialweakg'


nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    concured.append(float(item)) 
    #print(float(item))  
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[2] for data in nuevo0]:  
    cohered.append(float(item))  
    #print(float(item)) 
nuevo0.close()

#datos = 'parcialweakg'
datos1 = 'globalcouplarge'
#datos = 'parcialweakg'


nuevo1 = open(datos1)
for item in [data.split()[1] for data in nuevo1]:  
    concuglobal.append(float(item)) 
    #print(float(item))  
nuevo1.close()

nuevo1 = open(datos1)
for item in [data.split()[2] for data in nuevo1]:  
    coheglobal.append(float(item))  
    #print(float(item)) 
nuevo1.close()

plt.plot(g_data,coheparcial,label = "parcial", color = 'b',lw = 2)
plt.plot(g_data,cohered, linestyle='--', label = "redfield", color = 'r',lw=2)  
plt.plot(g_data,coheglobal, linestyle='--', label = "global", color = 'black',lw=2)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel("Coherencia",fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()


plt.plot(g_data,concuparcial,label = "parcial", color = 'b',lw = 2)
plt.plot(g_data,concured, linestyle='--', label = "redfield", color = 'r',lw=2)  
plt.plot(g_data,concuglobal, linestyle='--', label = "global", color = 'black',lw=2)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel("Concurrencia",fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()

