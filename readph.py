import os
import matplotlib.pyplot as plt

eV_data = []
Nparcial = []
Jlparcial = []
Jrparcial = []
JlRparcial = []
Jphparcial = []
coheparcial = []
concuparcial = []
Idparcial = []

Nred = []
Jlred = []
Jrred = []
JlRred = []
Jphred = []
cohered = []
concured = []
Idred = []

Nglobal = []
Jlglobal = []
Jrglobal = []
JlRglobal = []
Jphglobal = []
coheglobal = []
concuglobal = []
Idglobal = []
#######################################
########compararcodigos################
#datos = 'parcialphg'
datos = 'semiphg0'
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
    Jphparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()


nuevo = open(datos)
for item in [data.split()[6] for data in nuevo]:  
    coheparcial.append(float(item))   
nuevo.close()
nuevo = open(datos)
for item in [data.split()[7] for data in nuevo]:  
    concuparcial.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[8] for data in nuevo]:  
    Idparcial.append(float(item))   
nuevo.close()

#datos0 = 'redfieldweakg'
#datos0 = 'redfieldstrongg'
#datos0 = 'redfielphg'
datos0 = 'semiphJ0'
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
    Jphred.append(float(item))      
fichero0.close()

fichero0 = open(datos0)
for item in [data.split()[6] for data in fichero0]: #saca la columna 0|
    cohered.append(float(item))
fichero0.close()
fichero0 = open(datos0)
for item in [data.split()[7] for data in fichero0]: #saca la columna 0|
    concured.append(float(item))
fichero0.close()      
fichero0 = open(datos0)  
for item in [data.split()[8] for data in fichero0]: #saca la columna 0|
    Idred.append(float(item))
fichero0.close()  


#datosf = 'globalstrongg'
#datosf = 'globalphg'
datosf = 'semiphg0J0'

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
    Jphglobal.append(float(item))
ficherof.close()

ficherof = open(datosf)
for item in [data.split()[6] for data in ficherof]: #saca    
    coheglobal.append(float(item))
ficherof.close()

ficherof = open(datosf)
for item in [data.split()[7] for data in ficherof]: #saca    
    concuglobal.append(float(item))
ficherof.close()

ficherof = open(datosf)
for item in [data.split()[8] for data in ficherof]: #saca    
    Idglobal.append(float(item))
ficherof.close()


plt.plot(eV_data,Nparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,Nred, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,Nglobal, linestyle='--',color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jlparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,Jlred,color='blue',lw = 4, label = "J0")
plt.plot(eV_data,Jlglobal, linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Jrparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,Jrred, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,Jrglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{R}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,JlRparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,JlRred, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,JlRglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{LR}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,coheparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,cohered,color='blue',lw = 4,  label = "J0")
plt.plot(eV_data,coheglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$C$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,concuparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,concured, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,concuglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$Concu$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Idparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,Idred, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,Idglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{D}$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()


plt.plot(eV_data,Jphparcial, color='red',lw = 4, label = "g0")
plt.plot(eV_data,Jphred, color='blue',lw = 4, label = "J0")
plt.plot(eV_data,Jphglobal,linestyle='--', color='black',lw = 4, label = "g0,J0")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$J_{ph}$',fontsize = 25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()
