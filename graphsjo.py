import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phononJ0=10^{-3}bkb100.npz")
gof10,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0 = data0["gof1"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["concv"], data0["Nls"]

data1 = np.load("phononJ0=10^{-2}bkb100.npz")
gof11,Id1,Ile1,Ire1,Iphs1,cohes1,concv1,Nls1 = data1["gof1"], data1["Id"], data1["Ile"], data1["Ire"], data1["Iphs"], data1["cohes"], data1["concv"], data1["Nls"]

data2 = np.load("phononJ0=10^{-1}bkb100.npz")
gof12,Id2,Ile2,Ire2,Iphs2,cohes2,concv2,Nls2 = data2["gof1"], data2["Id"], data2["Ile"], data2["Ire"], data2["Iphs"], data2["cohes"], data2["concv"], data2["Nls"]

data3 = np.load("phononJ0=bkb100.npz")
gof13,Id3,Ile3,Ire3,Iphs3,cohes3,concv3,Nls3 = data3["gof1"], data3["Id"], data3["Ile"], data3["Ire"], data3["Iphs"], data3["cohes"], data3["concv"], data3["Nls"]

data4 = np.load("phononJ0=10bkb100.npz")
gof14,Id4,Ile4,Ire4,Iphs4,cohes4,concv4,Nls4 = data4["gof1"], data4["Id"], data4["Ile"], data4["Ire"], data4["Iphs"], data4["cohes"], data4["concv"], data4["Nls"]

data5 = np.load("phononJ0=5_10^{-1}bkb100.npz")
gof15,Id5,Ile5,Ire5,Iphs5,cohes5,concv5,Nls5 = data5["gof1"], data5["Id"], data5["Ile"], data5["Ire"], data5["Iphs"], data5["cohes"], data5["concv"], data5["Nls"]


Ilrt0 = []
Ilrt1 = []
Ilrt2 = []
Ilrt3 = []
Ilrt4 = []
Ilrt5 = []
Ilrt6 = []
Ilrt7 = []

for i in range(len(gof10)):
    Ilrt0.append(Ile0[i]+Ire0[i]+Iphs0[i])
    Ilrt1.append(Ile1[i]+Ire1[i]+Iphs1[i])
    Ilrt2.append(Ile2[i]+Ire2[i]+Iphs2[i])
    Ilrt3.append(Ile3[i]+Ire3[i]+Iphs3[i])
    Ilrt4.append(Ile4[i]+Ire4[i]+Iphs4[i])
    Ilrt5.append(Ile5[i]+Ire5[i]+Iphs5[i])
 

plt.plot(gof10,Nls0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}=10^{-3}$')
plt.plot(gof11,Nls1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}=10^{-2}$')
plt.plot(gof12,Nls2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,Nls3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,Nls4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,Nls5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(gof10,Id0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof11,Id1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof12,Id2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,Id3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,Id4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,Id5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{D}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(gof10,Ile0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof10,Ire0, color='blue',linestyle ='--',lw=3)
plt.plot(gof11,Ile1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof11,Ire1, color='orange',linestyle = '--',lw=3)
plt.plot(gof12,Ile2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof12,Ire2, color='green',linestyle = '--',lw=3)
plt.plot(gof13,Ile3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof13,Ire3, color='red',linestyle = '--',lw=3)
plt.plot(gof14,Ile4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof14,Ire4, color='purple',linestyle = '--',lw=3)
plt.plot(gof15,Ile5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.plot(gof15,Ire5, color='brown',linestyle = '--',lw=3)
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{le}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(gof10,concv0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof11,concv1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof12,concv2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,concv3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,concv4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,concv5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(gof10,cohes0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof11,cohes1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof12,cohes2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,cohes3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,cohes4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,cohes5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(gof10,Iphs0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof11,Iphs1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof12,Iphs2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,Iphs3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,Iphs4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,Iphs5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{ph}$',fontsize = 25) 
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(gof10,Ilrt0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
plt.plot(gof11,Ilrt1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
plt.plot(gof12,Ilrt2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
plt.plot(gof13,Ilrt3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
plt.plot(gof14,Ilrt4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
plt.plot(gof15,Ilrt5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
plt.xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{LRT}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns


#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax10.plot(gof10,Nls0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
ax10.plot(gof11,Nls1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
ax10.plot(gof12,Nls2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
ax10.plot(gof13,Nls3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
ax10.plot(gof14,Nls4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
ax10.plot(gof15,Nls5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
ax10.set_ylabel(r'$\dot{N}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "upper right")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

ax20.plot(gof10,concv0, color='blue',lw=3, label = r'$g=0$')
ax20.plot(gof11,concv1, color='orange',linestyle = '--',lw=3, label = r'$g=10^{-4}$')
ax20.plot(gof12,concv2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
ax20.plot(gof13,concv3, color='red',lw=3, label = r'$g=10^{-3}$')
ax20.plot(gof14,concv4, color='purple',lw=3,linestyle = '--', label = r'$g=5\cdot10^{-3}$')
ax20.plot(gof15,concv5, color='brown',lw=3, label = r'$g=5\cdot10^{-1}$')
ax20.set_xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Entrelazamiento",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()


# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns


#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax10.plot(gof10,Nls0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
ax10.plot(gof11,Nls1, color='orange',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
ax10.plot(gof12,Nls2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
ax10.plot(gof13,Nls3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
ax10.plot(gof14,Nls4, color='purple',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
ax10.plot(gof15,Nls5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
ax10.set_ylabel(r'$\dot{N}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "upper right")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

ax20.plot(gof10,cohes0, color='blue',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-3}$')
ax20.plot(gof11,cohes1, color='orange',linestyle = '--',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-2}$')
ax20.plot(gof12,cohes2, color='green',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10^{-1}$')
ax20.plot(gof13,cohes3, color='red',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 1$')
ax20.plot(gof14,cohes4, color='purple',lw=3,linestyle = '--', label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 10$')
ax20.plot(gof15,cohes5, color='brown',lw=3, label = r'$\frac{J0}{\beta_{ph} \kappa_{L}}= 5 \cdot 10^{-1}$')
ax20.set_xlabel(r'$g_{0}/(\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 

plt.plot(gof10,cohes0, color='blue',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 10^{-3}$')
plt.plot(gof11,cohes1, color='orange',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 10^{-2}$')
plt.plot(gof12,cohes2, color='green',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 10^{-1}$')
plt.plot(gof15,cohes5, color='brown',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 5 \cdot 10^{-1}$')
plt.plot(gof13,cohes3, color='red',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 1$')
plt.plot(gof14,cohes4, color='purple',lw=3, label = r'$J_{0}/\beta_{Ph} \kappa_{L}= 10$')
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 35)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 35)
plt.xticks(fontsize=35)  # X-axis tick labels
plt.yticks(fontsize=35)
plt.legend(fontsize=20,loc = "upper right")
plt.xscale("log")
plt.show()