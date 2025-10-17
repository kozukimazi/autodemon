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

Nglobal = []
Jlglobal = []
Jrglobal = []
JlRglobal = []
coheglobal = []
concuglobal = []

#######################################
########compararcodigos################
datos = 'parcialweakg'
#datos = 'parcialstrongg'
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
#datos0 = 'redfieldstrongg'
datos0 = 'redfieldweakg'
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

#datosf = 'globalstrongg'
datosf = 'globalweakg'

ficherof = open(datosf)
for item in [data.split()[1] for data in ficherof]: #saca       
    Nglobal.append(-float(item))
ficherof.close()
ficherof = open(datosf)
for item in [data.split()[2] for data in ficherof]: #saca       
    Jlglobal.append(float(item))
ficherof.close()

ficherof = open(datosf)
for item in [data.split()[3] for data in ficherof]: #saca    
    Jrglobal.append(float(item))
ficherof.close()
ficherof = open(datosf)
for item in [data.split()[4] for data in ficherof]: #saca    
    JlRglobal.append(float(item))
ficherof.close()
ficherof = open(datosf)
for item in [data.split()[5] for data in ficherof]: #saca    
    coheglobal.append(float(item))
ficherof.close()

ficherof = open(datosf)
for item in [data.split()[6] for data in ficherof]: #saca    
    concuglobal.append(float(item))
ficherof.close()



plt.plot(eV_data,Nparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Nred,linestyle='--', color='blue',lw = 4, label = "redfield")
plt.plot(eV_data,Nglobal, linestyle='--',color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jlparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Jlred, linestyle='--',color='blue',lw = 4, label = "redfield")
plt.plot(eV_data,Jlglobal, linestyle='--', color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jrparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,Jrred,linestyle='--', color='blue',lw = 4, label = "redfield")
plt.plot(eV_data,Jrglobal,linestyle='--', color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{R}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,JlRparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,JlRred,linestyle='--', color='blue',lw = 4, label = "redfield")
plt.plot(eV_data,JlRglobal,linestyle='--', color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{LR}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,coheparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,cohered, linestyle='--',color='blue',lw = 4,  label = "redfield")
plt.plot(eV_data,coheglobal,linestyle='--', color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$C$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,concuparcial, color='red',lw = 4, label = "parcial")
plt.plot(eV_data,concured,linestyle='--', color='blue',lw = 4, label = "redfield")
plt.plot(eV_data,concuglobal,linestyle='--', color='black',lw = 4, label = "global")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$Concu$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()
