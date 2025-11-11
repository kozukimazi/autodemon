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

plt.plot(Jof0,Nls0, color='blue',lw=3, label = r'$g=0$')
plt.plot(Jof1,Nls1, color='orange',lw=3, label = r'$g=10^{-4}$')
plt.plot(Jof2,Nls2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
plt.plot(Jof3,Nls3, color='red',lw=3, label = r'$g=10^{-3}$')
plt.plot(Jof4,Nls4, color='purple',lw=3, label = r'$g=5\cdot10^{-3}$')
plt.plot(Jof5,Nls5, color='brown',lw=3, label = r'$g=10^{-2}$')
plt.plot(Jof6,Nls6, color='pink',lw=3, label = r'$g=5\cdot10^{-2}$')
plt.plot(Jof7,Nls7, color='gray',lw=3, label = r'$g=10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Id0, color='blue',lw=3, label = r'$g=0$')
plt.plot(Jof1,Id1, color='orange',lw=3, label = r'$g=10^{-4}$')
plt.plot(Jof2,Id2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
plt.plot(Jof3,Id3, color='red',lw=3, label = r'$g=10^{-3}$')
plt.plot(Jof4,Id4, color='purple',lw=3, label = r'$g=5\cdot10^{-3}$')
plt.plot(Jof5,Id5, color='brown',lw=3, label = r'$g=10^{-2}$')
plt.plot(Jof6,Id6, color='pink',lw=3, label = r'$g=5\cdot10^{-2}$')
plt.plot(Jof7,Id7, color='gray',lw=3, label = r'$g=10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{D}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,Ile0, color='blue',lw=3, label = r'$g=0$')
plt.plot(Jof0,Ire0, color='blue',linestyle ='--',lw=3)
plt.plot(Jof1,Ile1, color='orange',lw=3, label = r'$g=10^{-4}$')
plt.plot(Jof1,Ire1, color='orange',linestyle = '--',lw=3)
plt.plot(Jof2,Ile2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
plt.plot(Jof2,Ire2, color='green',linestyle = '--',lw=3)
plt.plot(Jof3,Ile3, color='red',lw=3, label = r'$g=10^{-3}$')
plt.plot(Jof3,Ire3, color='red',linestyle = '--',lw=3)
plt.plot(Jof4,Ile4, color='purple',lw=3, label = r'$g=5\cdot10^{-3}$')
plt.plot(Jof4,Ire4, color='purple',linestyle = '--',lw=3)
plt.plot(Jof5,Ile5, color='brown',lw=3, label = r'$g=10^{-2}$')
plt.plot(Jof5,Ire5, color='brown',linestyle = '--',lw=3)
plt.plot(Jof6,Ile6, color='pink',lw=3, label = r'$g=5\cdot10^{-2}$')
plt.plot(Jof6,Ire6, color='pink',linestyle = '--',lw=3)
plt.plot(Jof7,Ile7, color='gray',lw=3, label = r'$g=10^{-1}$')
plt.plot(Jof7,Ire7, color='gray',linestyle = '--',lw=3)
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{le}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,concv0, color='blue',lw=3, label = r'$g=0$')
plt.plot(Jof1,concv1, color='orange',lw=3, label = r'$g=10^{-4}$')
plt.plot(Jof2,concv2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
plt.plot(Jof3,concv3, color='red',lw=3, label = r'$g=10^{-3}$')
plt.plot(Jof4,concv4, color='purple',lw=3, label = r'$g=5\cdot10^{-3}$')
plt.plot(Jof5,concv5, color='brown',lw=3, label = r'$g=10^{-2}$')
plt.plot(Jof6,concv6, color='pink',lw=3, label = r'$g=5\cdot10^{-2}$')
plt.plot(Jof7,concv7, color='gray',lw=3, label = r'$g=10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,cohes0, color='blue',lw=3, label = r'$g=0$')
plt.plot(Jof1,cohes1, color='orange',lw=3, label = r'$g=10^{-4}$')
plt.plot(Jof2,cohes2, color='green',lw=3, label = r'$g=5\cdot10^{-4}$')
plt.plot(Jof3,cohes3, color='red',lw=3, label = r'$g=10^{-3}$')
plt.plot(Jof4,cohes4, color='purple',lw=3, label = r'$g=5\cdot10^{-3}$')
plt.plot(Jof5,cohes5, color='brown',lw=3, label = r'$g=10^{-2}$')
plt.plot(Jof6,cohes6, color='pink',lw=3, label = r'$g=5\cdot10^{-2}$')
plt.plot(Jof7,cohes7, color='gray',lw=3, label = r'$g=10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()