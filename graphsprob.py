import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phonong=0prob.npz")
Jof0,Probnt10,Probnt20,Probnt30,Probnt40,Probnt50,Probnt60,Probnt70,Probnt80,Imalphg0,Imbetg0 = data0["Jof"], data0["Probnt10"], data0["Probnt20"], data0["Probnt30"], data0["Probnt40"], data0["Probnt50"], data0["Probnt60"], data0["Probnt70"], data0["Probnt80"], data0["Imalphg"], data0["Imbetg"]

data1 = np.load("phonong=10^{-4}prob.npz")
Jof1,Probnt11,Probnt21,Probnt31,Probnt41,Probnt51,Probnt61,Probnt71,Probnt81,Imalphg1,Imbetg1 = data1["Jof"], data1["Probnt10"], data1["Probnt20"], data1["Probnt30"], data1["Probnt40"], data1["Probnt50"], data1["Probnt60"], data1["Probnt70"], data1["Probnt80"], data1["Imalphg"], data1["Imbetg"]

data2 = np.load("phonong=5_10^{-4}prob.npz")
Jof2,Probnt12,Probnt22,Probnt32,Probnt42,Probnt52,Probnt62,Probnt72,Probnt82,Imalphg2,Imbetg2 = data2["Jof"], data2["Probnt10"], data2["Probnt20"], data2["Probnt30"], data2["Probnt40"], data2["Probnt50"], data2["Probnt60"], data2["Probnt70"], data2["Probnt80"], data2["Imalphg"], data2["Imbetg"]
data3 = np.load("phonong=10^{-3}prob.npz")
Jof3,Probnt13,Probnt23,Probnt33,Probnt43,Probnt53,Probnt63,Probnt73,Probnt83,Imalphg3,Imbetg3 = data3["Jof"], data3["Probnt10"], data3["Probnt20"], data3["Probnt30"], data3["Probnt40"], data3["Probnt50"], data3["Probnt60"], data3["Probnt70"], data3["Probnt80"], data3["Imalphg"], data3["Imbetg"]

data4 = np.load("phonong=5_10^{-3}prob.npz")
Jof4,Probnt14,Probnt24,Probnt34,Probnt44,Probnt54,Probnt64,Probnt74,Probnt84,Imalphg4,Imbetg4 = data4["Jof"], data4["Probnt10"], data4["Probnt20"], data4["Probnt30"], data4["Probnt40"], data4["Probnt50"], data4["Probnt60"], data4["Probnt70"], data4["Probnt80"], data4["Imalphg"], data4["Imbetg"]

data5 = np.load("phonong=10^{-2}prob.npz")
Jof5,Probnt15,Probnt25,Probnt35,Probnt45,Probnt55,Probnt65,Probnt75,Probnt85,Imalphg5,Imbetg5 = data5["Jof"], data5["Probnt10"], data5["Probnt20"], data5["Probnt30"], data5["Probnt40"], data5["Probnt50"], data5["Probnt60"], data5["Probnt70"], data5["Probnt80"], data5["Imalphg"], data5["Imbetg"]

data6 = np.load("phonong=5_10^{-2}prob.npz")
Jof6,Probnt16,Probnt26,Probnt36,Probnt46,Probnt56,Probnt66,Probnt76,Probnt86,Imalphg6,Imbetg6 = data6["Jof"], data6["Probnt10"], data6["Probnt20"], data6["Probnt30"], data6["Probnt40"], data6["Probnt50"], data6["Probnt60"], data6["Probnt70"], data6["Probnt80"], data6["Imalphg"], data6["Imbetg"]

data7 = np.load("phonong=10^{-1}prob.npz")
Jof7,Probnt17,Probnt27,Probnt37,Probnt47,Probnt57,Probnt67,Probnt77,Probnt87,Imalphg7,Imbetg7 = data7["Jof"], data7["Probnt10"], data7["Probnt20"], data7["Probnt30"], data7["Probnt40"], data7["Probnt50"], data7["Probnt60"], data7["Probnt70"], data7["Probnt80"], data7["Imalphg"], data7["Imbetg"]

data8 = np.load("phonong=7_10^{-3}prob.npz")
Jof8,Probnt18,Probnt28,Probnt38,Probnt48,Probnt58,Probnt68,Probnt78,Probnt88,Imalphg8,Imbetg8 = data8["Jof"], data8["Probnt10"], data8["Probnt20"], data8["Probnt30"], data8["Probnt40"], data8["Probnt50"], data8["Probnt60"], data8["Probnt70"], data8["Probnt80"], data8["Imalphg"], data8["Imbetg"]

data9 = np.load("phonong=3_10^{-3}prob.npz")
Jof9,Probnt19,Probnt29,Probnt39,Probnt49,Probnt59,Probnt69,Probnt79,Probnt89,Imalphg9,Imbetg9 = data9["Jof"], data9["Probnt10"], data9["Probnt20"], data9["Probnt30"], data9["Probnt40"], data9["Probnt50"], data9["Probnt60"], data9["Probnt70"], data9["Probnt80"], data9["Imalphg"], data9["Imbetg"]


plt.plot(Jof0,Probnt10, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt11, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt12, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt13, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt14, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt15, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt16, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt17, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof8,Probnt18, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Probnt19, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{111}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Probnt20, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt21, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt22, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt23, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt24, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt25, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt26, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt27, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt28, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Probnt29, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{110}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,Probnt30, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt31, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt32, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt33, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt34, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt35, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt36, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt37, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt38, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Probnt39, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{101}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Probnt40, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt41, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt42, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt43, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt44, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt45, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt46, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt47, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt48, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Probnt49, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{100}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Probnt50, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt51, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt52, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt53, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt54, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt55, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt56, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt57, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt58, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Probnt59, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{011}$',fontsize = 25) 
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Probnt60, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt61, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt62, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt63, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt64, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt65, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt66, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt67, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt68, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot10^{-1}$')
plt.plot(Jof9,Probnt69, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot10^{-1}$')

plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{010}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,Probnt70, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt71, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt72, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt73, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt74, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt75, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt76, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt77, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt78, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Probnt79, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{001}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Probnt80, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Probnt81, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Probnt82, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Probnt83, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Probnt84, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Probnt85, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Probnt86, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Probnt87, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Probnt88, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Probnt89, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$\rho_{000}$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof0,Imalphg0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Imalphg1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Imalphg2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Imalphg3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Imalphg4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Imalphg5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Imalphg6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Imalphg7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Imalphg8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Imalphg9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$2g\text{Im}(\alpha^{*})$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof0,Imbetg0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,Imbetg1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,Imbetg2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,Imbetg3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof4,Imbetg4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,Imbetg5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
plt.plot(Jof6,Imbetg6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,Imbetg7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
plt.plot(Jof8,Imbetg8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
plt.plot(Jof9,Imbetg9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)
plt.ylabel(r'$2g\text{Im}(\beta^{*})$',fontsize = 25)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()



