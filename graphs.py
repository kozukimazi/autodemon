import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phonong=0.npz")
Jof0,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0 = data0["Jof"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["concv"], data0["Nls"]

data1 = np.load("phonong=10^{-4}.npz")
Jof1,Id1,Ile1,Ire1,Iphs1,cohes1,concv1,Nls1 = data1["Jof"], data1["Id"], data1["Ile"], data1["Ire"], data1["Iphs"], data1["cohes"], data1["concv"], data1["Nls"]

data2 = np.load("phonong=5_10^{-4}.npz")
Jof2,Id2,Ile2,Ire2,Iphs2,cohes2,concv2,Nls2 = data2["Jof"], data2["Id"], data2["Ile"], data2["Ire"], data2["Iphs"], data2["cohes"], data2["concv"], data2["Nls"]

data3 = np.load("phonong=10^{-3}.npz")
Jof3,Id3,Ile3,Ire3,Iphs3,cohes3,concv3,Nls3 = data3["Jof"], data3["Id"], data3["Ile"], data3["Ire"], data3["Iphs"], data3["cohes"], data3["concv"], data3["Nls"]

data4 = np.load("phonong=5_10^{-3}.npz")
Jof4,Id4,Ile4,Ire4,Iphs4,cohes4,concv4,Nls4 = data4["Jof"], data4["Id"], data4["Ile"], data4["Ire"], data4["Iphs"], data4["cohes"], data4["concv"], data4["Nls"]

data5 = np.load("phonong=10^{-2}.npz")
Jof5,Id5,Ile5,Ire5,Iphs5,cohes5,concv5,Nls5 = data5["Jof"], data5["Id"], data5["Ile"], data5["Ire"], data5["Iphs"], data5["cohes"], data5["concv"], data5["Nls"]

data6 = np.load("phonong=5_10^{-2}.npz")
Jof6,Id6,Ile6,Ire6,Iphs6,cohes6,concv6,Nls6 = data6["Jof"], data6["Id"], data6["Ile"], data6["Ire"], data6["Iphs"], data6["cohes"], data6["concv"], data6["Nls"]

data7 = np.load("phonong=10^{-1}.npz")
Jof7,Id7,Ile7,Ire7,Iphs7,cohes7,concv7,Nls7 = data7["Jof"], data7["Id"], data7["Ile"], data7["Ire"], data7["Iphs"], data7["cohes"], data7["concv"], data7["Nls"]

data8 = np.load("phonong=7_10^{-3}.npz")
Jof8,Id8,Ile8,Ire8,Iphs8,cohes8,concv8,Nls8 = data8["Jof"], data8["Id"], data8["Ile"], data8["Ire"], data8["Iphs"], data8["cohes"], data8["concv"], data8["Nls"]

data9 = np.load("phonong=3_10^{-3}.npz")
Jof9,Id9,Ile9,Ire9,Iphs9,cohes9,concv9,Nls9 = data9["Jof"], data9["Id"], data9["Ile"], data9["Ire"], data9["Iphs"], data9["cohes"], data9["concv"], data9["Nls"]


Ilrt0 = []
Ilrt1 = []
Ilrt2 = []
Ilrt3 = []
Ilrt4 = []
Ilrt5 = []
Ilrt6 = []
Ilrt7 = []
Ilrt8 = []
Ilrt9 = []

for i in range(len(Jof0)):
    Ilrt0.append(Ile0[i]+Ire0[i])
    Ilrt1.append(Ile1[i]+Ire1[i])
    Ilrt2.append(Ile2[i]+Ire2[i])
    Ilrt3.append(Ile3[i]+Ire3[i])
    Ilrt4.append(Ile4[i]+Ire4[i])
    Ilrt5.append(Ile5[i]+Ire5[i])
    Ilrt6.append(Ile6[i]+Ire6[i])
    Ilrt7.append(Ile7[i]+Ire7[i])  
    Ilrt8.append(Ile8[i]+Ire8[i])  
    Ilrt9.append(Ile9[i]+Ire9[i])

plt.plot(Jof0,Nls0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Nls1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Nls2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Nls3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Nls4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Nls5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Nls6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Nls7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,Nls8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Nls9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Id0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Id1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Id2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Id3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Id4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Id5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Id6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Id7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Id8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Id9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{D}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,Ile0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof0,Ire0, color='blue',linestyle ='--',lw=3)
plt.plot(Jof1,Ile1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof1,Ire1, color='orange',linestyle = '--',lw=3)
plt.plot(Jof2,Ile2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof2,Ire2, color='green',linestyle = '--',lw=3)
plt.plot(Jof3,Ile3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof3,Ire3, color='red',linestyle = '--',lw=3)
plt.plot(Jof4,Ile4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof4,Ire4, color='purple',linestyle = '--',lw=3)
plt.plot(Jof5,Ile5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof5,Ire5, color='brown',linestyle = '--',lw=3)
plt.plot(Jof6,Ile6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof6,Ire6, color='pink',linestyle = '--',lw=3)
plt.plot(Jof7,Ile7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof7,Ire7, color='gray',linestyle = '--',lw=3)
plt.plot(Jof8,Ile8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof8,Ire8, color='black',linestyle = '--',lw=3)
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{le}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,concv0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,concv1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,concv2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,concv3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,concv4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,concv5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,concv6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,concv7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,concv8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,concv9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,cohes0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,cohes1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,cohes2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,cohes3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,cohes4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,cohes5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,cohes6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,cohes7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,cohes8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,cohes9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Iphs0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Iphs1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Iphs2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Iphs3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Iphs4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Iphs5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Iphs6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Iphs7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Iphs8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Iphs9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{ph}$',fontsize = 25) 
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Ilrt0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Ilrt1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Ilrt2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Ilrt3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Ilrt4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Ilrt5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Ilrt6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Ilrt7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Ilrt8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Ilrt9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
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
ax10.plot(Jof0,Nls0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax10.plot(Jof1,Nls1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax10.plot(Jof2,Nls2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax10.plot(Jof3,Nls3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax10.plot(Jof4,Nls4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax10.plot(Jof5,Nls5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax10.plot(Jof6,Nls6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax10.plot(Jof7,Nls7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
ax10.plot(Jof8,Nls8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
ax10.plot(Jof9,Nls9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

ax10.set_ylabel(r'$\dot{N}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "upper right")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

ax20.plot(Jof0,concv0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax20.plot(Jof1,concv1, color='orange',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax20.plot(Jof2,concv2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax20.plot(Jof3,concv3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax20.plot(Jof4,concv4, color='purple',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax20.plot(Jof5,concv5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax20.plot(Jof6,concv6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax20.plot(Jof7,concv7, color='gray',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=10^{1}$') 
ax20.plot(Jof8,concv8, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$') 
ax20.plot(Jof9,concv9, color='black',lw=3,linestyle = ':', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$') 

ax20.set_xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia y entrelazamiento",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()


# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns


#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax10.plot(Jof0,Nls0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax10.plot(Jof1,Nls1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax10.plot(Jof2,Nls2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax10.plot(Jof3,Nls3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax10.plot(Jof4,Nls4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax10.plot(Jof5,Nls5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax10.plot(Jof6,Nls6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax10.plot(Jof7,Nls7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
ax10.plot(Jof8,Nls8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
ax10.plot(Jof9,Nls9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

ax10.set_ylabel(r'$\dot{N}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "upper right")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

ax20.plot(Jof0,cohes0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax20.plot(Jof1,cohes1, color='orange',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax20.plot(Jof2,cohes2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax20.plot(Jof3,cohes3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax20.plot(Jof4,cohes4, color='purple',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax20.plot(Jof5,cohes5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax20.plot(Jof6,cohes6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax20.plot(Jof7,cohes7, color='gray',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=10^{1}$') 
ax20.plot(Jof8,cohes8, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$') 
ax20.plot(Jof9,cohes9, color='black',lw=3,linestyle = ':', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$') 

ax20.set_xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()

# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns

ax10.plot(Jof0,Iphs0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax10.plot(Jof1,Iphs1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax10.plot(Jof2,Iphs2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax10.plot(Jof3,Iphs3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax10.plot(Jof4,Iphs4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax10.plot(Jof5,Iphs5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax10.plot(Jof6,Iphs6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax10.plot(Jof7,Iphs7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
ax10.plot(Jof8,Iphs8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
ax10.plot(Jof9,Iphs9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
ax10.set_ylabel(r'$\dot{I}_{ph}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "upper right")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


ax20.plot(Jof0,cohes0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
ax20.plot(Jof1,cohes1, color='orange',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
ax20.plot(Jof2,cohes2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
ax20.plot(Jof3,cohes3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
ax20.plot(Jof4,cohes4, color='purple',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
ax20.plot(Jof5,cohes5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
ax20.plot(Jof6,cohes6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
ax20.plot(Jof7,cohes7, color='gray',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=10^{1}$') 
ax20.plot(Jof8,cohes8, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$') 
ax20.plot(Jof9,cohes9, color='black',lw=3,linestyle = ':', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$') 
ax20.set_xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()
