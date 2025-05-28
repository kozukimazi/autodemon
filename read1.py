import os
import matplotlib.pyplot as plt


eV_data = []
Jclas_data = []
Prob1clas = []
Prob2clas = []
Prob3clas = []
Prob4clas = []
Prob5clas = []
Prob6clas = []
Prob7clas = []
Prob8clas = []
Iclas_data = []


Jqm_data = []
Prob1qm = []
Prob2qm = []
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
    Prob1clas.append(float(item))  
    #print(float(item)) 
nuevo.close()

nuevo = open(datos)
for item in [data.split()[3] for data in nuevo]:  
    Prob2clas.append(float(item))  
    #print(float(item)) 
nuevo.close()
        
nuevo = open(datos)
for item in [data.split()[4] for data in nuevo]:  
    Prob3clas.append(float(item))  
    #print(float(item)) 
nuevo.close()
nuevo = open(datos)
for item in [data.split()[5] for data in nuevo]:  
    Prob4clas.append(float(item))   
nuevo.close()
nuevo = open(datos)
for item in [data.split()[6] for data in nuevo]:  
    Prob5clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[7] for data in nuevo]:  
    Prob6clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[8] for data in nuevo]:  
    Prob7clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[9] for data in nuevo]:  
    Prob8clas.append(float(item))   
nuevo.close()

nuevo = open(datos)
for item in [data.split()[10] for data in nuevo]:  
    Iclas_data.append(float(item))   
nuevo.close()


datos0 = 'lindbladgamU'
nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    Jqm_data.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)     
for item in [data.split()[2] for data in nuevo0]:  
    Prob1qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)     
for item in [data.split()[3] for data in nuevo0]:  
    Prob2qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)     
for item in [data.split()[4] for data in nuevo0]:  
    Prob3qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[5] for data in nuevo0]:  
    Prob4qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[6] for data in nuevo0]:  
    Prob5qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[7] for data in nuevo0]:  
    Prob6qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[8] for data in nuevo0]:  
    Prob7qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[9] for data in nuevo0]:  
    Prob8qm.append(float(item))   
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[10] for data in nuevo0]:  
    Iqm_data.append(float(item))   
nuevo0.close()

plt.plot(eV_data,Jclas_data, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Jqm_data, color='blue',lw = 4, label = "quantum")
plt.plot(eV_data,aux0)
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob1clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob1qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{101}$',fontsize = 25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15)
plt.show()

plt.plot(eV_data,Prob2clas, color='red',lw = 4, label = "classic")
plt.plot(eV_data,Prob2qm, color='blue',lw = 4, label = "quantum")
plt.xlabel(r'$eV/T$',fontsize = 20)
plt.ylabel(r'$\rho_{101}$',fontsize = 25)
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


# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# First plot
ax1.plot(eV_data,Prob4clas, color='red',lw = 4, label = "classic")
ax1.plot(eV_data,Prob4qm, color='blue',lw = 4, label = "quantum")
ax1.set_xlabel(r'$eV/T$',fontsize=20)
ax1.set_ylabel(r'$\rho_{100}$',fontsize = 20)
ax1.legend(fontsize=15)
ax1.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Second plot
ax2.plot(eV_data,Prob6clas, color='red',lw = 4, label = "classic")
ax2.plot(eV_data,Prob6qm, color='blue',lw = 4, label = "quantum")
ax2.set_xlabel(r'$eV/T$',fontsize=20)
ax2.set_ylabel(r'$\rho_{010}$',fontsize = 20)
ax2.legend(fontsize=15)
ax2.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plots
plt.show()


# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# First plot
ax1.plot(eV_data,Prob3clas, color='red',lw = 4, label = "classic")
ax1.plot(eV_data,Prob3qm, color='blue',lw = 4, label = "quantum")
ax1.set_xlabel(r'$eV/T$',fontsize=20)
ax1.set_ylabel(r'$\rho_{101}$',fontsize = 20)
ax1.legend(fontsize=13.8,loc = 'upper left')
ax1.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Second plot
ax2.plot(eV_data,Prob5clas, color='red',lw = 4, label = "classic")
ax2.plot(eV_data,Prob5qm, color='blue',lw = 4, label = "quantum")
ax2.set_xlabel(r'$eV/T$',fontsize=20)
ax2.set_ylabel(r'$\rho_{011}$',fontsize = 20)
ax2.legend(fontsize=14,loc = 'upper left')
ax2.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plots
plt.show()

# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# First plot
ax1.plot(eV_data,Prob7clas, color='red',lw = 4, label = "classic")
ax1.plot(eV_data,Prob7qm, color='blue',lw = 4, label = "quantum")
ax1.set_xlabel(r'$eV/T$',fontsize=20)
ax1.set_ylabel(r'$\rho_{001}$',fontsize = 20)
ax1.legend(fontsize=14,loc = 'upper right')
ax1.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Second plot
ax2.plot(eV_data,Prob8clas, color='red',lw = 4, label = "classic")
ax2.plot(eV_data,Prob8qm, color='blue',lw = 4, label = "quantum")
ax2.set_xlabel(r'$eV/T$',fontsize=20)
ax2.set_ylabel(r'$\rho_{000}$',fontsize = 20)
ax2.legend(fontsize=14,loc = 'upper right')
ax2.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plots
plt.show()


# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# First plot
ax1.plot(eV_data,Jclas_data, color='red',lw = 4, label = "classic")
ax1.plot(eV_data,Jqm_data, color='blue',lw = 4, label = "quantum")
ax1.set_xlabel(r'$eV/T$',fontsize=20)
ax1.set_ylabel(r'$\dot{N}_{L}$',fontsize = 20)
ax1.legend(fontsize=14,loc = 'upper left')
ax1.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Second plot
ax2.plot(eV_data,Iclas_data, color='red',lw = 4, label = "classic")
ax2.plot(eV_data,Iqm_data, color='blue',lw = 4, label = "quantum")
ax2.set_xlabel(r'$eV/T$',fontsize=20)
ax2.set_ylabel(r'$\dot{I}_{D}$',fontsize = 20)
ax2.legend(fontsize=14,loc = 'upper right')
ax2.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plots
plt.show()


# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# First plot
ax1.plot(eV_data,Prob1clas, color='red',lw = 4, label = "classic")
ax1.plot(eV_data,Prob1qm, color='blue',lw = 4, label = "quantum")
ax1.set_xlabel(r'$eV/T$',fontsize=20)
ax1.set_ylabel(r'$\rho_{111}$',fontsize = 20)
ax1.legend(fontsize=14,loc = 'upper left')
ax1.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Second plot
ax2.plot(eV_data,Prob2clas, color='red',lw = 4, label = "classic")
ax2.plot(eV_data,Prob2qm, color='blue',lw = 4, label = "quantum")
ax2.set_xlabel(r'$eV/T$',fontsize=20)
ax2.set_ylabel(r'$\rho_{110}$',fontsize = 20)
ax2.legend(fontsize=14,loc = 'upper left')
ax2.tick_params(axis='both', which='major', labelsize=10)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=8)   # Minor ticks

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plots
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