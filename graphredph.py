import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phonong=0red.npz")
Jof0,cohev0,concuv0,Nls0,Jlrs0,Wlr0,Jphs0,Imalph0,Imbetg0 = data0["Jof"], data0["cohev"], data0["concuv"], data0["Nls"], data0["Jlrs"], data0["Wlr"], data0["Jphs"], data0["Imalphg"], data0["Imbetg"]

data1 = np.load("phonong=10^{-4}red.npz")
Jof1,cohev1,concuv1,Nls1,Jlrs1,Wlr1,Jphs1,Imalph1,Imbetg1 = data1["Jof"], data1["cohev"], data1["concuv"], data1["Nls"], data1["Jlrs"], data1["Wlr"], data1["Jphs"], data1["Imalphg"], data1["Imbetg"]

data2 = np.load("phonong=5_10^{-4}red.npz")
Jof2,cohev2,concuv2,Nls2,Jlrs2,Wlr2,Jphs2,Imalph2,Imbetg2 = data2["Jof"], data2["cohev"], data2["concuv"], data2["Nls"], data2["Jlrs"], data2["Wlr"], data2["Jphs"], data2["Imalphg"], data2["Imbetg"]
data3 = np.load("phonong=10^{-3}red.npz")
Jof3,cohev3,concuv3,Nls3,Jlrs3,Wlr3,Jphs3,Imalph3,Imbetg3 = data3["Jof"], data3["cohev"], data3["concuv"], data3["Nls"], data3["Jlrs"], data3["Wlr"], data3["Jphs"], data3["Imalphg"], data3["Imbetg"]

data4 = np.load("phonong=5_10^{-3}red.npz")
Jof4,cohev4,concuv4,Nls4,Jlrs4,Wlr4,Jphs4,Imalph4,Imbetg4 = data4["Jof"], data4["cohev"], data4["concuv"], data4["Nls"], data4["Jlrs"], data4["Wlr"], data4["Jphs"], data4["Imalphg"], data4["Imbetg"]

data5 = np.load("phonong=10^{-2}red.npz")
Jof5,cohev5,concuv5,Nls5,Jlrs5,Wlr5,Jphs5,Imalph5,Imbetg5 = data5["Jof"], data5["cohev"], data5["concuv"], data5["Nls"], data5["Jlrs"], data5["Wlr"], data5["Jphs"], data5["Imalphg"], data5["Imbetg"]

data6 = np.load("phonong=5_10^{-2}red.npz")
Jof6,cohev6,concuv6,Nls6,Jlrs6,Wlr6,Jphs6,Imalph6,Imbetg6 = data6["Jof"], data6["cohev"], data6["concuv"], data6["Nls"], data6["Jlrs"], data6["Wlr"], data6["Jphs"], data6["Imalphg"], data6["Imbetg"]

data7 = np.load("phonong=10^{-1}red.npz")
Jof7,cohev7,concuv7,Nls7,Jlrs7,Wlr7,Jphs7,Imalph7,Imbetg7 = data7["Jof"], data7["cohev"], data7["concuv"], data7["Nls"], data7["Jlrs"], data7["Wlr"], data7["Jphs"], data7["Imalphg"], data7["Imbetg"]

data8 = np.load("phonong=7_10^{-3}red.npz")
Jof8,cohev8,concuv8,Nls8,Jlrs8,Wlr8,Jphs8,Imalph8,Imbetg8 = data8["Jof"], data8["cohev"], data8["concuv"], data8["Nls"], data8["Jlrs"], data8["Wlr"], data8["Jphs"], data8["Imalphg"], data8["Imbetg"]

data9 = np.load("phonong=3_10^{-3}red.npz")
Jof9,cohev9,concuv9,Nls9,Jlrs9,Wlr9,Jphs9,Imalph9,Imbetg9 = data9["Jof"], data9["cohev"], data9["concuv"], data9["Nls"], data9["Jlrs"], data9["Wlr"], data9["Jphs"], data9["Imalphg"], data9["Imbetg"]

data10 = np.load("phonong=1red.npz")
Jof10,cohev10,concuv10,Nls10,Jlrs10,Wlr10,Jphs10,Imalph10,Imbetg10 = data10["Jof"], data10["cohev"], data10["concuv"], data10["Nls"], data10["Jlrs"], data10["Wlr"], data10["Jphs"], data10["Imalphg"], data10["Imbetg"]



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
plt.plot(Jof10,Nls10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{N}_{L}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,cohev0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,cohev1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,cohev2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,cohev3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,cohev4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,cohev5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,cohev6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,cohev7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,cohev8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,cohev9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.plot(Jof10,cohev10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,concuv0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,concuv1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,concuv2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,concuv3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,concuv4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,concuv5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,concuv6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,concuv7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,concuv8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,concuv9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.plot(Jof10,concuv10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Jlrs0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Jlrs1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Jlrs2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Jlrs3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Jlrs4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Jlrs5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Jlrs6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Jlrs7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,Jlrs8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Jlrs9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.plot(Jof10,Jlrs10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$J_{LRS}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Wlr0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Wlr1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Wlr2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Wlr3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Wlr4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Wlr5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Wlr6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Wlr7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,Wlr8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Wlr9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.plot(Jof10,Wlr10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$W_{LR}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Imalph0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Imalph1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Imalph2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Imalph3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Imalph4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Imalph5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Imalph6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Imalph7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,Imalph8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Imalph9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.plot(Jof10,Imalph10, color='cyan',lw=3, label = r'$\frac{g}{\kappa_{L}}=100$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$2g \text{Im}(\alpha)$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

