import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext,NullLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


data0 = np.load("phonong=0probb100.npz")
Jof0,Probnt10,Probnt20,Probnt30,Probnt40,Probnt50,Probnt60,Probnt70,Probnt80,Imalphg0,Imbetg0 = data0["Jof"], data0["Probnt10"], data0["Probnt20"], data0["Probnt30"], data0["Probnt40"], data0["Probnt50"], data0["Probnt60"], data0["Probnt70"], data0["Probnt80"], data0["Imalphg"], data0["Imbetg"]



data1 = np.load("phonong=10^{-4}probb100.npz")
Jof1,Probnt11,Probnt21,Probnt31,Probnt41,Probnt51,Probnt61,Probnt71,Probnt81,Imalphg1,Imbetg1 = data1["Jof"], data1["Probnt10"], data1["Probnt20"], data1["Probnt30"], data1["Probnt40"], data1["Probnt50"], data1["Probnt60"], data1["Probnt70"], data1["Probnt80"], data1["Imalphg"], data1["Imbetg"]

data2 = np.load("phonong=5_10^{-4}probb100.npz")
Jof2,Probnt12,Probnt22,Probnt32,Probnt42,Probnt52,Probnt62,Probnt72,Probnt82,Imalphg2,Imbetg2 = data2["Jof"], data2["Probnt10"], data2["Probnt20"], data2["Probnt30"], data2["Probnt40"], data2["Probnt50"], data2["Probnt60"], data2["Probnt70"], data2["Probnt80"], data2["Imalphg"], data2["Imbetg"]
data3 = np.load("phonong=10^{-3}probb100.npz")
Jof3,Probnt13,Probnt23,Probnt33,Probnt43,Probnt53,Probnt63,Probnt73,Probnt83,Imalphg3,Imbetg3 = data3["Jof"], data3["Probnt10"], data3["Probnt20"], data3["Probnt30"], data3["Probnt40"], data3["Probnt50"], data3["Probnt60"], data3["Probnt70"], data3["Probnt80"], data3["Imalphg"], data3["Imbetg"]

data4 = np.load("phonong=5_10^{-3}probb100.npz")
Jof4,Probnt14,Probnt24,Probnt34,Probnt44,Probnt54,Probnt64,Probnt74,Probnt84,Imalphg4,Imbetg4 = data4["Jof"], data4["Probnt10"], data4["Probnt20"], data4["Probnt30"], data4["Probnt40"], data4["Probnt50"], data4["Probnt60"], data4["Probnt70"], data4["Probnt80"], data4["Imalphg"], data4["Imbetg"]

data5 = np.load("phonong=10^{-2}probb100.npz")
Jof5,Probnt15,Probnt25,Probnt35,Probnt45,Probnt55,Probnt65,Probnt75,Probnt85,Imalphg5,Imbetg5 = data5["Jof"], data5["Probnt10"], data5["Probnt20"], data5["Probnt30"], data5["Probnt40"], data5["Probnt50"], data5["Probnt60"], data5["Probnt70"], data5["Probnt80"], data5["Imalphg"], data5["Imbetg"]

data6 = np.load("phonong=5_10^{-2}probb100.npz")
Jof6,Probnt16,Probnt26,Probnt36,Probnt46,Probnt56,Probnt66,Probnt76,Probnt86,Imalphg6,Imbetg6 = data6["Jof"], data6["Probnt10"], data6["Probnt20"], data6["Probnt30"], data6["Probnt40"], data6["Probnt50"], data6["Probnt60"], data6["Probnt70"], data6["Probnt80"], data6["Imalphg"], data6["Imbetg"]

data7 = np.load("phonong=10^{-1}probb100.npz")
Jof7,Probnt17,Probnt27,Probnt37,Probnt47,Probnt57,Probnt67,Probnt77,Probnt87,Imalphg7,Imbetg7 = data7["Jof"], data7["Probnt10"], data7["Probnt20"], data7["Probnt30"], data7["Probnt40"], data7["Probnt50"], data7["Probnt60"], data7["Probnt70"], data7["Probnt80"], data7["Imalphg"], data7["Imbetg"]

data8 = np.load("phonong=7_10^{-3}probb100.npz")
Jof8,Probnt18,Probnt28,Probnt38,Probnt48,Probnt58,Probnt68,Probnt78,Probnt88,Imalphg8,Imbetg8 = data8["Jof"], data8["Probnt10"], data8["Probnt20"], data8["Probnt30"], data8["Probnt40"], data8["Probnt50"], data8["Probnt60"], data8["Probnt70"], data8["Probnt80"], data8["Imalphg"], data8["Imbetg"]

data9 = np.load("phonong=3_10^{-3}probb100.npz")
Jof9,Probnt19,Probnt29,Probnt39,Probnt49,Probnt59,Probnt69,Probnt79,Probnt89,Imalphg9,Imbetg9 = data9["Jof"], data9["Probnt10"], data9["Probnt20"], data9["Probnt30"], data9["Probnt40"], data9["Probnt50"], data9["Probnt60"], data9["Probnt70"], data9["Probnt80"], data9["Imalphg"], data9["Imbetg"]

data10 = np.load("phonong=10probb100.npz")
Jof10,Probnt110,Probnt210,Probnt310,Probnt410,Probnt510,Probnt610,Probnt710,Probnt810,Imalphg10,Imbetg10 = data10["Jof"], data10["Probnt10"], data10["Probnt20"], data10["Probnt30"], data10["Probnt40"], data10["Probnt50"], data10["Probnt60"], data10["Probnt70"], data10["Probnt80"], data10["Imalphg"], data10["Imbetg"]

data11 = np.load("phonong=100probb100.npz")
Jof11,Probnt111,Probnt211,Probnt311,Probnt411,Probnt511,Probnt611,Probnt711,Probnt811,Imalphg11,Imbetg11 = data11["Jof"], data11["Probnt10"], data11["Probnt20"], data11["Probnt30"], data11["Probnt40"], data11["Probnt50"], data11["Probnt60"], data11["Probnt70"], data11["Probnt80"], data11["Imalphg"], data11["Imbetg"]

data12 = np.load("phonong=1000probb100.npz")
Jof12,Probnt112,Probnt212,Probnt312,Probnt412,Probnt512,Probnt612,Probnt712,Probnt812,Imalphg12,Imbetg12 = data12["Jof"], data12["Probnt10"], data12["Probnt20"], data12["Probnt30"], data12["Probnt40"], data12["Probnt50"], data12["Probnt60"], data12["Probnt70"], data12["Probnt80"], data12["Imalphg"], data12["Imbetg"]

data13 = np.load("phonong=10000probb100.npz")
Jof13,Probnt113,Probnt213,Probnt313,Probnt413,Probnt513,Probnt613,Probnt713,Probnt813,Imalphg13,Imbetg13 = data13["Jof"], data13["Probnt10"], data13["Probnt20"], data13["Probnt30"], data13["Probnt40"], data13["Probnt50"], data13["Probnt60"], data13["Probnt70"], data13["Probnt80"], data13["Imalphg"], data13["Imbetg"]

data14 = np.load("phonong=1probb100.npz")
Jof14,Probnt114,Probnt214,Probnt314,Probnt414,Probnt514,Probnt614,Probnt714,Probnt814,Imalphg14,Imbetg14 = data14["Jof"], data14["Probnt10"], data14["Probnt20"], data14["Probnt30"], data14["Probnt40"], data14["Probnt50"], data14["Probnt60"], data14["Probnt70"], data14["Probnt80"], data14["Imalphg"], data14["Imbetg"]

dataf0 = np.load("phonong=0b100.npz")
Jof0f,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0,Act0,Nqm0,Nltotal0,Work0,eff0,effph0,coheveig0 = dataf0["Jof"], dataf0["Id"], dataf0["Ile"], dataf0["Ire"], dataf0["Iphs"], dataf0["cohes"], dataf0["concv"], dataf0["Nls"],dataf0["Acts"], dataf0["Nlqm"],dataf0["Nltotal"],dataf0["Work"],dataf0["eff"],dataf0["effph"],dataf0["coheveig"]

dataf1 = np.load("phonong=10^{-4}b100.npz")
Jof1f,Id1,Ile1,Ire1,Iphs1,cohes1,concv1,Nls1,Act1,Nqm1,Nltotal1,Work1,eff1,effph1,coheveig1 = dataf1["Jof"], dataf1["Id"], dataf1["Ile"], dataf1["Ire"], dataf1["Iphs"], dataf1["cohes"], dataf1["concv"], dataf1["Nls"],dataf1["Acts"], dataf1["Nlqm"],dataf1["Nltotal"],dataf1["Work"],dataf1["eff"],dataf1["effph"],dataf1["coheveig"]   

dataf2 = np.load("phonong=5_10^{-4}b100.npz")
Jof2f,Id2,Ile2,Ire2,Iphs2,cohes2,concv2,Nls2,Act2,Nqm2,Nltotal2,Work2,eff2,effph2,coheveig2 = dataf2["Jof"], dataf2["Id"], dataf2["Ile"], dataf2["Ire"], dataf2["Iphs"], dataf2["cohes"], dataf2["concv"], dataf2["Nls"],dataf2["Acts"], dataf2["Nlqm"],dataf2["Nltotal"],dataf2["Work"],dataf2["eff"],dataf2["effph"],dataf2["coheveig"]

dataf3 = np.load("phonong=10^{-3}b100.npz")
Jof3f,Id3,Ile3,Ire3,Iphs3,cohes3,concv3,Nls3,Act3,Nqm3,Nltotal3,Work3,eff3,effph3,coheveig3 = dataf3["Jof"], dataf3["Id"], dataf3["Ile"], dataf3["Ire"], dataf3["Iphs"], dataf3["cohes"], dataf3["concv"], dataf3["Nls"],dataf3["Acts"], dataf3["Nlqm"],dataf3["Nltotal"],dataf3["Work"],dataf3["eff"],dataf3["effph"],dataf3["coheveig"]

dataf4 = np.load("phonong=5_10^{-3}b100.npz")
Jof4f,Id4,Ile4,Ire4,Iphs4,cohes4,concv4,Nls4,Act4,Nqm4,Nltotal4,Work4,eff4,effph4,coheveig4 = dataf4["Jof"], dataf4["Id"], dataf4["Ile"], dataf4["Ire"], dataf4["Iphs"], dataf4["cohes"], dataf4["concv"], dataf4["Nls"],dataf4["Acts"], dataf4["Nlqm"],dataf4["Nltotal"],dataf4["Work"],dataf4["eff"],dataf4["effph"],dataf4["coheveig"]

dataf6 = np.load("phonong=5_10^{-2}b100.npz")
Jof6f,Id6,Ile6,Ire6,Iphs6,cohes6,concv6,Nls6,Act6,Nqm6,Nltotal6,Work6,eff6,effph6,coheveig6 = dataf6["Jof"], dataf6["Id"], dataf6["Ile"], dataf6["Ire"], dataf6["Iphs"], dataf6["cohes"], dataf6["concv"], dataf6["Nls"],dataf6["Acts"], dataf6["Nlqm"],dataf6["Nltotal"],dataf6["Work"],dataf6["eff"],dataf6["effph"],dataf6["coheveig"]


dataf7 = np.load("phonong=10^{-1}b100.npz")
Jof7f,Id7,Ile7,Ire7,Iphs7,cohes7,concv7,Nls7,Act7,Nqm7,Nltotal7,Work7,eff7,effph7,coheveig7 = dataf7["Jof"], dataf7["Id"], dataf7["Ile"], dataf7["Ire"], dataf7["Iphs"], dataf7["cohes"], dataf7["concv"], dataf7["Nls"],dataf7["Acts"], dataf7["Nlqm"],dataf7["Nltotal"],dataf7["Work"],dataf7["eff"],dataf7["effph"],dataf7["coheveig"]


dataf8 = np.load("phonong=7_10^{-3}b100.npz")
Jof8f,Id8,Ile8,Ire8,Iphs8,cohes8,concv8,Nls8,Act8,Nqm8,Nltotal8,Work8,eff8,effph8,coheveig8 = dataf8["Jof"], dataf8["Id"], dataf8["Ile"], dataf8["Ire"], dataf8["Iphs"], dataf8["cohes"], dataf8["concv"], dataf8["Nls"],dataf8["Acts"], dataf8["Nlqm"],dataf8["Nltotal"],dataf8["Work"],dataf8["eff"],dataf8["effph"],dataf8["coheveig"]
dataf9 = np.load("phonong=3_10^{-3}b100.npz")
Jof9f,Id9,Ile9,Ire9,Iphs9,cohes9,concv9,Nls9,Act9,Nqm9,Nltotal9,Work9,eff9,effph9,coheveig9 = dataf9["Jof"], dataf9["Id"], dataf9["Ile"], dataf9["Ire"], dataf9["Iphs"], dataf9["cohes"], dataf9["concv"], dataf9["Nls"],dataf9["Acts"], dataf9["Nlqm"],dataf9["Nltotal"],dataf9["Work"],dataf9["eff"],dataf9["effph"],dataf9["coheveig"]

dataf10 = np.load("phonong=10b100.npz")
Jof10f,Id10,Ile10,Ire10,Iphs10,cohes10,concv10,Nls10,Work10,eff10,effph10,coheveig10 = dataf10["Jof"], dataf10["Id"], dataf10["Ile"], dataf10["Ire"], dataf10["Iphs"], dataf10["cohes"], dataf10["concv"], dataf10["Nls"],dataf10["Work"],dataf10["eff"],dataf10["effph"],dataf10["coheveig"]    

dataf11 = np.load("phonong=100b100.npz")
Jof11f,Id11,Ile11,Ire11,Iphs11,cohes11,concv11,Nls11,Work11,eff11,effph11,coheveig11 = dataf11["Jof"], dataf11["Id"], dataf11["Ile"], dataf11["Ire"], dataf11["Iphs"], dataf11["cohes"], dataf11["concv"], dataf11["Nls"],dataf11["Work"],dataf11["eff"],dataf11["effph"],dataf11["coheveig"]

dataf12 = np.load("phonong=1000b100.npz")
Jof12f,Id12,Ile12,Ire12,Iphs12,cohes12,concv12,Nls12,Work12,eff12,effph12,coheveig12 = dataf12["Jof"], dataf12["Id"], dataf12["Ile"], dataf12["Ire"], dataf12["Iphs"], dataf12["cohes"], dataf12["concv"], dataf12["Nls"],dataf12["Work"],dataf12["eff"],dataf12["effph"],dataf12["coheveig"]

dataf13 = np.load("phonong=10000b100.npz")
Jof13f,Id13,Ile13,Ire13,Iphs13,cohes13,concv13,Nls13,Work13,eff13,effph13,coheveig13 = dataf13["Jof"], dataf13["Id"], dataf13["Ile"], dataf13["Ire"], dataf13["Iphs"], dataf13["cohes"], dataf13["concv"], dataf13["Nls"],dataf13["Work"],dataf13["eff"],dataf13["effph"],dataf13["coheveig"]

dataf14 = np.load("phonong=1b100.npz")  
Jof14f,Id14,Ile14,Ire14,Iphs14,cohes14,concv14,Nls14,Work14,eff14,effph14,coheveig14 = dataf14["Jof"], dataf14["Id"], dataf14["Ile"], dataf14["Ire"], dataf14["Iphs"], dataf14["cohes"], dataf14["concv"], dataf14["Nls"],dataf14["Work"],dataf14["eff"],dataf14["effph"],dataf14["coheveig"]

dataf15 = np.load("phonong=5b100.npz")
Jof15f,Id15,Ile15,Ire15,Iphs15,cohes15,concv15,Nls15,Work15,eff15,effph15,coheveig15 = dataf15["Jof"], dataf15["Id"], dataf15["Ile"], dataf15["Ire"], dataf15["Iphs"], dataf15["cohes"], dataf15["concv"], dataf15["Nls"],dataf15["Work"],dataf15["eff"],dataf15["effph"],dataf15["coheveig"]

dataf16 = np.load("phonong=9b100.npz")
Jof16f,Id16,Ile16,Ire16,Iphs16,cohes16,concv16,Nls16,Work16,eff16,effph16,coheveig16 = dataf16["Jof"], dataf16["Id"], dataf16["Ile"], dataf16["Ire"], dataf16["Iphs"], dataf16["cohes"], dataf16["concv"], dataf16["Nls"],dataf16["Work"],dataf16["eff"],dataf16["effph"],dataf16["coheveig"]
######################


dataf5 = np.load("phonong=10^{-2}b100.npz")
Jof5f,Id5,Ile5,Ire5,Iphs5,cohes5,concv5,Nls5,Act5,Nqm5,Nltotal5,Work5,eff5,effph5,coheveig5 = dataf5["Jof"], dataf5["Id"], dataf5["Ile"], dataf5["Ire"], dataf5["Iphs"], dataf5["cohes"], dataf5["concv"], dataf5["Nls"],dataf5["Acts"], dataf5["Nlqm"],dataf5["Nltotal"],dataf5["Work"],dataf5["eff"],dataf5["effph"],dataf5["coheveig"]


n = len(Imalphg0)
cohesum0s = []
cohesum1s = []
cohesum2s = []
cohesum3s = []
cohesum4s = []
cohesum5s = []
cohesum6s = []
cohesum7s = []
cohesum8s = []
cohesum9s = []
cohesum10s = []
cohesum11s = []
cohesum12s = []
cohesum13s = []
cohesum14s = []
for i in range(n):
    Nls0[i] = abs(Nls0[i]*100)
    Nls1[i] = abs(Nls1[i]*100)
    Nls2[i] = abs(Nls2[i]*100)
    Nls3[i] = abs(Nls3[i]*100)
    Nls4[i] = abs(Nls4[i]*100)
    Nls5[i] = abs(Nls5[i]*100)
    Nls6[i] = abs(Nls6[i]*100)
    Nls7[i] = abs(Nls7[i]*100)
    Nls8[i] = abs(Nls8[i]*100)
    Nls9[i] = abs(Nls9[i]*100)
    Nls10[i] = abs(Nls10[i]*100)
    Nls11[i] = abs(Nls11[i]*100)
    Nls12[i] = abs(Nls12[i]*100)
    Nls13[i] = abs(Nls13[i]*100)
    Nls14[i] = abs(Nls14[i]*100)
    Nls15[i] = abs(Nls15[i]*100)
    cohesum0s.append((Imalphg0[i] + Imbetg0[i])*100)
    cohesum1s.append((Imalphg1[i] + Imbetg1[i])*100)
    cohesum2s.append((Imalphg2[i] + Imbetg2[i])*100)
    cohesum3s.append((Imalphg3[i] + Imbetg3[i])*100)
    cohesum4s.append((Imalphg4[i] + Imbetg4[i])*100)
    cohesum5s.append((Imalphg5[i] + Imbetg5[i])*100)
    cohesum6s.append((Imalphg6[i] + Imbetg6[i])*100)    
    cohesum7s.append((Imalphg7[i] + Imbetg7[i])*100)
    cohesum8s.append((Imalphg8[i] + Imbetg8[i])*100)
    cohesum9s.append((Imalphg9[i] + Imbetg9[i])*100)
    cohesum10s.append((Imalphg10[i] + Imbetg10[i])*100)
    cohesum11s.append((Imalphg11[i] + Imbetg11[i])*100)
    cohesum12s.append((Imalphg12[i] + Imbetg12[i])*100)
    cohesum13s.append((Imalphg13[i] + Imbetg13[i])*100)
    cohesum14s.append((Imalphg14[i] + Imbetg14[i])*100)

'''

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

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 

'''


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 


# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(3.39, 4.8))  # 1 row, 2 columns
TICK_FS = 8

#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax20.plot(Jof0,Imalphg0, color='blue',lw=1.6, label = r'$g/\kappa_{L}=0$')
#ax20.plot(Jof1,Imalphg1, color='orange',lw=1.6, label = r'$g/\kappa_{L}=10^{-2}$')
ax20.plot(Jof2,Imalphg2, color='green',lw=1.6, label = r'$g/\kappa_{L}=5\times 10^{-2}$')
ax20.plot(Jof3,Imalphg3, color='red',lw=1.6, label = r'$g/\kappa_{L}=10^{-1}$')
ax20.plot(Jof9,Imalphg9, color='black',lw=1.6, label = r'$g/\kappa_{L}=3\times 10^{-1}$')
ax20.plot(Jof4,Imalphg4, color='purple',lw=1.6, label = r'$g/\kappa_{L}=5\times 10^{-1}$')
ax20.plot(Jof5,Imalphg5, color='brown',lw=1.6, label = r'$g/\kappa_{L}=10^{0}$')
#plt.plot(Jof6,Imalphg6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
#ax20.plot(Jof7,Imalphg7, color='gray',lw=1.6, label = r'$g/\kappa_{L}=10^{1}$')
#plt.plot(Jof8,Imalphg8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')

ax20.set_ylabel(r'$2g|\mathrm{Im}(\alpha)|$',fontsize = 9)
ax20.set_xlabel(r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_{L})$',fontsize = 9) 
ax20.set_xscale("log")  
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.9, 0.97, '(b)', transform=ax20.transAxes, fontsize=9, fontweight='bold', va='top', ha='left')

ax10.plot(Jof0,Imbetg0, color='blue',lw=1.6, label = r'$g/\kappa_{L}=0$')
#ax10.plot(Jof1,Imbetg1, color='orange',lw=1.6, label = r'$g/\kappa_{L}=10^{-2}$')
ax10.plot(Jof2,Imbetg2, color='green',lw=1.6, label = r'$g/\kappa_{L}=5\times 10^{-2}$')
ax10.plot(Jof3,Imbetg3, color='red',lw=1.6, label = r'$g/\kappa_{L}=10^{-1}$')
ax10.plot(Jof9,Imbetg9, color='black',lw=1.6, label = r'$g/\kappa_{L}=3 \times 10^{-1}$')
ax10.plot(Jof4,Imbetg4, color='purple',lw=1.6, label = r'$g/\kappa_{L}=5 \times 10^{-1}$')
ax10.plot(Jof5,Imbetg5, color='brown',lw=1.6, label = r'$g/\kappa_{L}=10^{0}$')
#plt.plot(Jof6,Imbetg6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
#ax10.plot(Jof7,Imbetg7, color='gray',lw=1.6, label = r'$g/\kappa_{L}=10^{1}$')
#plt.plot(Jof8,Imbetg8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
ax10.legend(fontsize=8,loc = "center left", bbox_to_anchor=(0.02, 0.65))  
  
ax10.set_xscale("log")
#ax20.legend(fontsize=17, loc = "upper left") 
ax10.set_ylabel(r'$2g|\mathrm{Im}(\beta)|$',fontsize = 9)
ax10.tick_params(direction='in', which='both',labelbottom = False, labelsize=TICK_FS)  # font size of tick labels 
ax10.text(0.9, 0.97, '(a)', transform=ax10.transAxes, fontsize=9, fontweight='bold', va='top', ha='left')

ax10.xaxis.set_major_locator(LogLocator(base=10))
ax10.xaxis.set_minor_locator(LogLocator(base=10, subs=[2,5]))
ax10.xaxis.set_major_formatter(LogFormatterMathtext())
# Remove ALL minor ticks
ax10.xaxis.set_minor_locator(NullLocator())

ax20.xaxis.set_major_locator(LogLocator(base=10))
ax20.xaxis.set_minor_locator(LogLocator(base=10, subs=[2,5]))
ax20.xaxis.set_major_formatter(LogFormatterMathtext())
ax20.xaxis.set_minor_locator(NullLocator())

plt.tight_layout(pad=0.4)  # Avoids overlapping labels
plt.savefig("fig6_PRapp.pdf")
plt.show()


fig, (ax10, ax20) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(3.39, 4.0)
)

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# ---------- Panel (a)
ax10.plot(Jof0, Imbetg0, color='blue',   lw=LINE_W, label=r'$g/\kappa_L=0$')
#ax10.plot(Jof1, Imalphg1, color='orange', lw=LINE_W, label=r'$10^{-2}$')
ax10.plot(Jof2, Imbetg2, color='green',  lw=LINE_W, label=r'$g/\kappa_L=5\times10^{-2}$')
ax10.plot(Jof3, Imbetg3, color='red',    lw=LINE_W, label=r'$g/\kappa_L=10^{-1}$')
ax10.plot(Jof9, Imbetg9, color='black',  lw=LINE_W, ls='--', label=r'$g/\kappa_L=3\times10^{-1}$')
ax10.plot(Jof4, Imbetg4, color='purple', lw=LINE_W, label=r'$g/\kappa_L=5\times10^{-1}$')
#ax10.plot(Jof5, Imalphg5, color='brown',  lw=LINE_W, label=r'$1$')
ax10.plot(Jof7, Imbetg7, color='gray',   lw=LINE_W, label=r'$g/\kappa_L=10$')
ax10.plot(Jof12, Imbetg12, color='green',   lw=LINE_W,ls='--', label=r'$g/\kappa_L=10^{5}$')
ax10.plot(Jof13, Imbetg13, color='brown',   lw=LINE_W,ls = '--', label=r'$g/\kappa_L=10^{6}$')

ax10.set_ylabel(r'$2g\mathrm{Im}(\beta)$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='center left'
)

# ---------- Panel (b)
ax20.plot(Jof0, Imalphg0, color='blue',   lw=LINE_W)
#ax20.plot(Jof1, Imbetg1, color='orange', lw=LINE_W)
ax20.plot(Jof2, Imalphg2, color='green',  lw=LINE_W)
ax20.plot(Jof3, Imalphg3, color='red',    lw=LINE_W)
ax20.plot(Jof9, Imalphg9, color='black',  lw=LINE_W, ls='--')
ax20.plot(Jof4, Imalphg4, color='purple', lw=LINE_W)
#ax20.plot(Jof5, Imbetg5, color='brown',  lw=LINE_W)
ax20.plot(Jof7, Imalphg7, color='gray',   lw=LINE_W)
ax20.plot(Jof12, Imalphg12, color='green',   lw=LINE_W,ls='--')
ax20.plot(Jof13, Imalphg13, color='brown',   lw=LINE_W,ls = '--')

ax20.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}}\kappa_L)$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$2g\,\mathrm{Im}(\alpha)$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("fig6_PR_colors0.pdf")
plt.close()


datos0 = np.load("phonong=0.1probb100sec.npz")
Jofs0,Probnts10,Probnts20,Probnts30,Probnts40,Probnts50,Probnts60,Probnts70,Probnts80,Imalphg0s,Imbetg0s = datos0["Jof"], datos0["Probnt10"], datos0["Probnt20"], datos0["Probnt30"], datos0["Probnt40"], data0["Probnt50"], datos0["Probnt60"], datos0["Probnt70"], datos0["Probnt80"], datos0["Imalphg"], datos0["Imbetg"]

datos1 = np.load("phonong=1probb100sec.npz")
Jofs1,Probnts11,Probnts21,Probnts31,Probnts41,Probnts51,Probnts61,Probnts71,Probnts81,Imalphg1s,Imbetg1s = datos1["Jof"], datos1["Probnt10"], datos1["Probnt20"], datos1["Probnt30"], datos1["Probnt40"], datos1["Probnt50"], datos1["Probnt60"], datos1["Probnt70"], datos1["Probnt80"], datos1["Imalphg"], datos1["Imbetg"]

datos2 = np.load("phonong=10probb100sec.npz")
Jofs2,Probnts12,Probnts22,Probnts32,Probnts42,Probnts52,Probnts62,Probnts72,Probnts82,Imalphg2s,Imbetg2s = datos2["Jof"], datos2["Probnt10"], datos2["Probnt20"], datos2["Probnt30"], datos2["Probnt40"], datos2["Probnt50"], datos2["Probnt60"], datos2["Probnt70"], datos2["Probnt80"], datos2["Imalphg"], datos2["Imbetg"]

datos3 = np.load("phonong=100probb100sec.npz")
Jofs3,Probnts13,Probnts23,Probnts33,Probnts43,Probnts53,Probnts63,Probnts73,Probnts83,Imalphg3s,Imbetg3s = datos3["Jof"], datos3["Probnt10"], datos3["Probnt20"], datos3["Probnt30"], datos3["Probnt40"], datos3["Probnt50"], datos3["Probnt60"], datos3["Probnt70"], datos3["Probnt80"], datos3["Imalphg"], datos3["Imbetg"]

datos4 = np.load("phonong=1000probb100sec.npz")
Jofs4,Probnts14,Probnts24,Probnts34,Probnts44,Probnts54,Probnts64,Probnts74,Probnts84,Imalphg4s,Imbetg4s = datos4["Jof"], datos4["Probnt10"], datos4["Probnt20"], datos4["Probnt30"], datos4["Probnt40"], datos4["Probnt50"], datos4["Probnt60"], datos4["Probnt70"], datos4["Probnt80"], datos4["Imalphg"], datos4["Imbetg"]





fig, (ax10, ax20) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(3.39, 4.0)
)

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# ---------- Panel (a)
ax10.plot(Jofs0, Imbetg0s, color='blue',   lw=LINE_W, label=r'$g/\kappa_L=10$')
ax10.plot(Jofs1, Imbetg1s, color='orange', lw=LINE_W, label=r'$g/\kappa_L=10^{2}$')
ax10.plot(Jofs2, Imbetg2s, color='green',  lw=LINE_W, label=r'$g/\kappa_L=10^{3}$')
ax10.plot(Jofs3, Imbetg3s, color='red',    lw=LINE_W, label=r'$g/\kappa_L=10^{4}$')
ax10.plot(Jofs4, Imbetg4s, color='purple', lw=LINE_W, label=r'$g/\kappa_L=10^{5}$')
#ax10.plot(Jof5, Imalphg5, color='brown',  lw=LINE_W, label=r'$1$')

ax10.set_ylabel(r'$2g\mathrm{Re}(\beta)$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax10.text(0.9, 0.87, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='center left'
)

# ---------- Panel (b)
ax20.plot(Jofs0, Imalphg0s, color='blue',   lw=LINE_W)
ax20.plot(Jofs1, Imalphg1s, color='orange', lw=LINE_W)
ax20.plot(Jofs2, Imalphg2s, color='green',  lw=LINE_W)
ax20.plot(Jofs3, Imalphg3s, color='red',    lw=LINE_W)
ax20.plot(Jofs4, Imalphg4s, color='purple', lw=LINE_W)
#ax20.plot(Jof5, Imbetg5, color='brown',  lw=LINE_W)

ax20.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}}\kappa_L)$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$2g\,\mathrm{Re}(\alpha)$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.9, 0.87, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("fig6_PR_colors0sec.pdf")
plt.show()
plt.close()


fig, ax10 = plt.subplots(
    sharex=True,
    figsize=(3.39, 2.8)
)

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# ---------- Panel (a)
ax10.plot(Jof0, cohesum0s, color='blue',   lw=LINE_W, label=r'$g/\kappa_L=0$')
#ax10.plot(Jof1, cohesum1s, color='orange', lw=LINE_W, label=r'$g/\kappa_L=10^{-2}$')
ax10.plot(Jof2, cohesum2s, color='green',  lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{-2}$')
ax10.plot(Jof3, cohesum3s, color='red',    lw=LINE_W, label=r'$g/\kappa_L=10^{-1}$')
ax10.plot(Jof4, cohesum4s, color='purple',  lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{-1}$')
ax10.plot(Jof5, cohesum5s, color='brown', lw=LINE_W, label=r'$g/\kappa_L=10^{0}$')
#ax10.plot(Jof6, cohesum6s, color='pink',   lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{0}$')
#ax10.plot(Jof7, cohesum7s, color='gray',   lw=LINE_W, label=r'$g/\kappa_L=10^{1}$')
#ax10.plot(Jof8, cohesum8s, color='black',  lw=LINE_W, label   =r'$g/\kappa_L=7 \times 10^{-1}$')
ax10.plot(Jof9, cohesum9s, color='black',  lw=LINE_W, label=r'$g/\kappa_L=3 \times 10^{-1}$')
#ax10.plot(Jof5, Imalphg5, color='brown',  lw=LINE_W, label=r'$1$')
ax10.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}}\kappa_L)$', fontsize=LABEL_FS)
ax10.set_ylabel(r'$2g|\mathrm{Im}(\alpha)+\mathrm{Im}(\beta)|$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
#ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes,
          #fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=6.1,
    frameon=True,
    ncol=1,
    loc='center left',
    bbox_to_anchor=(0.02, 0.63)
)



# ---------- Spines

for spine in ax10.spines.values():
    spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("fig6cohesum.pdf")
plt.close()




fig, ax10 = plt.subplots(
    sharex=True,
    figsize=(3.39, 2.8)
)

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# ---------- Panel (a)
ax10.plot(Jof0, cohesum0s, color='blue',   lw=LINE_W, label=r'$g/\kappa_L=0$')
ax10.plot(Jof0,Nls0, color='blue',   lw=LINE_W, ls='--')
#ax10.plot(Jof1, cohesum1s, color='orange', lw=LINE_W, label=r'$g/\kappa_L=10^{-2}$')
ax10.plot(Jof2, cohesum2s, color='green',  lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{-2}$')
ax10.plot(Jof2,Nls2, color='green',  lw=LINE_W, ls='--')
ax10.plot(Jof3, cohesum3s, color='red',    lw=LINE_W, label=r'$g/\kappa_L=10^{-1}$')
ax10.plot(Jof3,Nls3, color='red',    lw=LINE_W, ls='--')
ax10.plot(Jof4, cohesum4s, color='purple',  lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{-1}$')
ax10.plot(Jof4,Nls4, color='purple',  lw=LINE_W, ls='--')
ax10.plot(Jof5, cohesum5s, color='brown', lw=LINE_W, label=r'$g/\kappa_L=10^{0}$')
ax10.plot(Jof5,Nls5, color='brown', lw=LINE_W, ls='--')
#ax10.plot(Jof6, cohesum6s, color='pink',   lw=LINE_W, label=r'$g/\kappa_L=5 \times 10^{0}$')
#ax10.plot(Jof7, cohesum7s, color='gray',   lw=LINE_W, label=r'$g/\kappa_L=10^{1}$')
#ax10.plot(Jof8, cohesum8s, color='black',  lw=LINE_W, label   =r'$g/\kappa_L=7 \times 10^{-1}$')
ax10.plot(Jof9, cohesum9s, color='black',  lw=LINE_W, label=r'$g/\kappa_L=3 \times 10^{-1}$')
ax10.plot(Jof9,Nls9, color='black',  lw=LINE_W, ls='--')
#ax10.plot(Jof5, Imalphg5, color='brown',  lw=LINE_W, label=r'$1$')
ax10.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}}\kappa_L)$', fontsize=LABEL_FS)
#ax10.set_ylabel(r'$2g|\mathrm{Im}(\alpha)+\mathrm{Im}(\beta)|$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
#ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes,
          #fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=6.1,
    frameon=True,
    ncol=1,
    loc='center left',
    bbox_to_anchor=(0.02, 0.63)
)



# ---------- Spines

for spine in ax10.spines.values():
    spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("fig6coheNls.pdf")
plt.close()



#########hacer los calculos de 1800 a 3000
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    sharex=True,
    figsize=(7.0, 3.6)   # APS wide figure
)

LINE_W = 2.0
LABEL_FS = 9
TICK_FS  = 8
PANEL_FS = 9
LEG_FS   = 7

# ===================== (a) Efficiency =====================
ax1.plot(Jof0f, Nls0, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=0$')
#ax1.plot(Jof1, Nls1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
ax1.plot(Jof2f, Nls2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\times 10^{-2}$')
ax1.plot(Jof3f, Nls3, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{-1}$')
ax1.plot(Jof9f, Nls9, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=3\times10^{-1}$')
ax1.plot(Jof4f, Nls4, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=5\times10^{-1}$')
ax1.plot(Jof5f, Nls5, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{0}$')
#ax1.plot(Jof7, Nls7, color='gray',   lw=LINE_W, label=r'$g/\kappa_{L}=10^{1}$')
#ax1.plot(Jof10, Nls10, color='red',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{3}$')
#ax1.plot(Jof11, Nls11, color='blue',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{4}$')
#ax1.plot(Jof12, Nls12, color='green',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{5}$')
#ax1.plot(Jof13, Nls13, color='brown',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{6}$')
ax1.plot(Jof14f, Nls14, color='black',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{2}$')
#ax1.plot(Jof15, Nls15, color='gray',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=5\times 10^{2}$')
#ax1.plot(Jof16, Nls16, color='pink',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=9\times 10^{2}$')

axins = inset_axes(
    ax1,
    width="50%",
    height="50%",
    bbox_to_anchor=(0.35, 0.05, 0.6, 0.6),
    bbox_transform=ax1.transAxes,
    
    borderpad=0
)

# Plot SAME curves (no legend!)
axins.plot(Jof4, Nls4, color='purple', lw=2)
axins.plot(Jof5,Nls5,color='brown', lw = 2)
axins.plot(Jof9, Nls9, color='black',  lw=2)

# Zoom region (example values — adjust!)
axins.set_xlim(6.5e-4, 2.1e-0)
axins.set_ylim(0.0385, 0.0485)

# Inset formatting
axins.tick_params(direction='in', labelsize=6)
for spine in axins.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)


ax1.set_xscale("log")
ax1.set_ylabel(r'$|\dot{N}_{B_L}|/\kappa_{L}$', fontsize=LABEL_FS)

ax1.text(0.92, 0.90, '(a)', transform=ax1.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
ax2.plot(Jof0, cohesum0s, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, coheveig1, color='orange', lw=LINE_W)
ax2.plot(Jof2, cohesum2s, color='green',  lw=LINE_W)
ax2.plot(Jof3, cohesum3s, color='red',    lw=LINE_W)
ax2.plot(Jof9, cohesum9s, color='black',  lw=LINE_W)
ax2.plot(Jof4, cohesum4s, color='purple', lw=LINE_W)
ax2.plot(Jof5, cohesum5s, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, cohesum7s, color='gray',   lw=LINE_W)
#ax2.plot(Jof10, cohesum10s, color='red',   lw=LINE_W, ls = '--')
#ax2.plot(Jof11, cohesum11s, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, cohesum12s, color='green',   lw=LINE_W, ls = '--')
#ax2.plot(Jof13f, cohesum13s, color='brown',   lw=LINE_W, ls = '--')
ax2.plot(Jof14f, cohesum14s, color='black',   lw=LINE_W, ls = '--')
#ax2.plot(Jof15, cohesum15s, color='gray',   lw=LINE_W, ls = '--')
#ax2.plot(Jof16, cohesum16s, color='pink',   lw=LINE_W, ls = '--')



ax2.set_xscale("log")
ax2.set_ylabel(r'$\langle \dot{N}_{\mathrm{cohe}}\rangle/\kappa_{L}$', fontsize=LABEL_FS)

ax2.text(0.92, 0.90, '(b)', transform=ax2.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== Shared formatting =====================
for ax in (ax1, ax2):
    ax.tick_params(direction='in', which='both', labelsize=TICK_FS)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

ax1.set_xlabel(r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_{L})$', fontsize=LABEL_FS)
ax2.set_xlabel(r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_{L})$', fontsize=LABEL_FS)

# ===================== Single legend (recommended) =====================
handles, labels = ax1.get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc='upper center',
    ncol=4,
    fontsize=LEG_FS,
    frameon=True
)

# ===================== Layout & save =====================
plt.tight_layout(pad=0.4, rect=[0, 0, 1, 0.88])
plt.savefig("figcurrentdecohe.pdf")
plt.close()



#########hacer los calculos de 1800 a 3000
fig, (ax10, ax20) = plt.subplots(
    1, 2,
    sharex=True,
    figsize=(7.0, 3.6)   # APS wide figure
)

LINE_W = 2.0
LABEL_FS = 9
TICK_FS  = 8
PANEL_FS = 9
LEG_FS   = 7

# ===================== (a) Efficiency =====================
ax10.plot(Jof0, Id0, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=0$')
#ax1.plot(Jof1, eff1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
ax10.plot(Jof2, Id2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\times 10^{-2}$')
ax10.plot(Jof3, Id3, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{-1}$')
ax10.plot(Jof9, Id9, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=3\times10^{-1}$')
ax10.plot(Jof4, Id4, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=5\times10^{-1}$')
ax10.plot(Jof5, Id5, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{0}$')
#ax1.plot(Jof7, eff7, color='gray',   lw=LINE_W, label=r'$g/\kappa_{L}=10^{1}$')
#ax10.plot(Jof10, Id10, color='red',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{3}$')
#ax10.plot(Jof11, Id11, color='blue',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{4}$')
#ax1.plot(Jof12, eff12, color='green',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{5}$')
#ax1.plot(Jof13, eff13, color='brown',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{6}$')
ax10.plot(Jof14, Id14, color='black',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{2}$')
#ax1.plot(Jof15, eff15, color='gray',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=5\times 10^{2}$')
#ax1.plot(Jof16, eff16, color='pink',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=9\times 10^{2}$')

axins0 = inset_axes(
    ax10,
    width="30%",
    height="30%",
    #bbox_to_anchor=(0.005, 0.005, 0.1, 0.2),  #[x, y, width, height] 
    #bbox_transform=ax10.transAxes,
    loc = 'center right',
    borderpad=1
)

# Plot SAME curves (no legend!)
axins0.plot(Jof4, Id4, color='purple', lw=2)
axins0.plot(Jof9, Id9, color='black',  lw=2)
axins0.plot(Jof5, Id5, color='brown',  lw=2)
axins0.set_xscale("log")

# Zoom region (example values — adjust!)
axins0.set_xlim(9e-3, 2)
axins0.set_ylim(0.0042, 0.00431)

# Inset formatting
axins0.tick_params(direction='in', labelsize=6)
for spine in axins0.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
#mark_inset(ax10, axins0, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)
mark_inset(ax10, axins0, loc1=2, loc2=4, fc="none", ec="0.4", lw=0.6)

ax10.set_xscale("log")
ax10.set_ylabel(r'$\dot{I}_{1}$', fontsize=LABEL_FS)

ax10.text(0.92, 0.90, '(a)', transform=ax10.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
ax20.plot(Jof0, eff0, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, coheveig1, color='orange', lw=LINE_W)
ax20.plot(Jof2, eff2, color='green',  lw=LINE_W)
ax20.plot(Jof3, eff3, color='red',    lw=LINE_W)
ax20.plot(Jof9, eff9, color='black',  lw=LINE_W)
ax20.plot(Jof4, eff4, color='purple', lw=LINE_W)
ax20.plot(Jof5, eff5, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, coheveig7, color='gray',   lw=LINE_W)
#ax20.plot(Jof10, eff10, color='red',   lw=LINE_W, ls = '--')
#ax20.plot(Jof11, eff11, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, coheveig12, color='green',   lw=LINE_W, ls = '--')
#ax20.plot(Jof13, eff13, color='brown',   lw=LINE_W, ls = '--')
ax20.plot(Jof14, eff14, color='black',   lw=LINE_W, ls = '--')
#ax2.plot(Jof15, coheveig15, color='gray',   lw=LINE_W, ls = '--')
#ax2.plot(Jof16, coheveig16, color='pink',   lw=LINE_W, ls = '--')



ax20.set_xscale("log")
ax20.set_ylabel(r'$\eta_{2}$', fontsize=LABEL_FS)

ax20.text(0.92, 0.90, '(b)', transform=ax20.transAxes,
         fontsize=PANEL_FS, fontweight='bold')


axins20 = inset_axes(
    ax20,
    width="32%",
    height="32%",
    #bbox_to_anchor=(0.005, 0.005, 0.1, 0.2),  #[x, y, width, height] 
    #bbox_transform=ax10.transAxes,
    loc = 'center right',
    borderpad=1
)

# Plot SAME curves (no legend!)
axins20.plot(Jof4, eff4, color='purple', lw=2)
axins20.plot(Jof9, eff9, color='black',  lw=2)
axins20.plot(Jof5, eff5, color='brown',  lw=2)
axins20.set_xscale("log")

# Zoom region (example values — adjust!)
axins20.set_xlim(1e-2, 4)
axins20.set_ylim(0.5, 0.541)

# Inset formatting
axins20.tick_params(direction='in', labelsize=6)
for spine in axins20.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
#mark_inset(ax10, axins0, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)
mark_inset(ax20, axins20, loc1=2, loc2=4, fc="none", ec="0.4", lw=0.6)



# ===================== Shared formatting =====================
for ax in (ax10, ax20):
    ax.tick_params(direction='in', which='both', labelsize=TICK_FS)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

ax10.set_xlabel(r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_{L})$', fontsize=LABEL_FS)
ax20.set_xlabel(r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_{L})$', fontsize=LABEL_FS)

# ===================== Single legend (recommended) =====================
handles, labels = ax10.get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc='upper center',
    ncol=4,
    fontsize=LEG_FS,
    frameon=True
)

# ===================== Layout & save =====================
plt.tight_layout(pad=0.4, rect=[0, 0, 1, 0.88])
#plt.savefig("figinfocoheeig.pdf")
plt.close()




import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# ==========================================================
# PRB STYLE
# ==========================================================
LINE_W  = 2.0
LABEL_FS = 9
TICK_FS  = 8
PANEL_FS = 10
LEG_FS   = 7

# ==========================================================
# CREATE 2x2 FIGURE
# ==========================================================
fig, axs = plt.subplots(
    2, 2,
    figsize=(7.0, 6.0),
    sharex=True
)

(ax1, ax2), (ax3, ax4) = axs

# ==========================================================
# COMMON CURVES
# ==========================================================
curves = [
    (Jof0,  'blue',   r'$g/\kappa_L=0$'),
    (Jof2,  'green',  r'$g/\kappa_L=5\times10^{-2}$'),
    (Jof3,  'red',    r'$g/\kappa_L=10^{-1}$'),
    (Jof9,  'black',  r'$g/\kappa_L=3\times10^{-1}$'),
    (Jof4,  'purple', r'$g/\kappa_L=5\times10^{-1}$'),
    (Jof5,  'brown',  r'$g/\kappa_L=10^{0}$'),
    (Jof14, 'black',  r'$g/\kappa_L=10^{2}$')
]

# ==========================================================
# (a) CURRENT
# ==========================================================
datasets_a = [Nls0, Nls2, Nls3, Nls9, Nls4, Nls5, Nls14]

for (x, c, lab), y in zip(curves, datasets_a):
    ls = '--' if '10^{2}' in lab else '-'
    ax1.plot(x, y, color=c, lw=LINE_W, ls=ls, label=lab)

ax1.set_xscale('log')
ax1.set_ylabel(r'$|\dot{N}_{B_L}|/\kappa_L$', fontsize=LABEL_FS)

ax1.text(
    0.90, 0.89, '(a)',
    transform=ax1.transAxes,
    fontsize=PANEL_FS,
    fontweight='bold'
)

# ==========================================================
# INSET (a)
# ==========================================================
axins1 = inset_axes(
    ax1,
    width="34%",
    height="34%",
    loc='center right',
    borderpad=1
)

axins1.plot(Jof4, Nls4, color='purple', lw=1.8)
axins1.plot(Jof5, Nls5, color='brown', lw=1.8)
axins1.plot(Jof9, Nls9, color='black', lw=1.8)

axins1.set_xscale("log")

axins1.set_xlim(6.5e-4, 2.1)
axins1.set_ylim(0.0385, 0.0485)

axins1.tick_params(direction='in', labelsize=5)

for spine in axins1.spines.values():
    spine.set_linewidth(0.6)

mark_inset(
    ax1,
    axins1,
    loc1=2,
    loc2=4,
    fc="none",
    ec="0.4",
    lw=0.6
)

# ==========================================================
# (b) COHERENCE
# ==========================================================
datasets_b = [
    cohesum0s,
    cohesum2s,
    cohesum3s,
    cohesum9s,
    cohesum4s,
    cohesum5s,
    cohesum14s
]

for (x, c, lab), y in zip(curves, datasets_b):
    ls = '--' if '10^{2}' in lab else '-'
    ax2.plot(x, y, color=c, lw=LINE_W, ls=ls)

ax2.set_xscale('log')

ax2.set_ylabel(
    r'$\langle \dot{N}_{\mathrm{cohe}}\rangle/\kappa_L$',
    fontsize=LABEL_FS
)

ax2.text(
    0.90, 0.89, '(b)',
    transform=ax2.transAxes,
    fontsize=PANEL_FS,
    fontweight='bold'
)

# ==========================================================
# (c) INFORMATION FLOW
# ==========================================================
datasets_c = [Id0, Id2, Id3, Id9, Id4, Id5, Id14]

for (x, c, lab), y in zip(curves, datasets_c):
    ls = '--' if '10^{2}' in lab else '-'
    ax3.plot(x, y, color=c, lw=LINE_W, ls=ls)

ax3.set_xscale('log')

ax3.set_ylabel(
    r'$\dot{I}_1$',
    fontsize=LABEL_FS
)

ax3.text(
    0.90, 0.89, '(c)',
    transform=ax3.transAxes,
    fontsize=PANEL_FS,
    fontweight='bold'
)

# ==========================================================
# INSET (c)
# ==========================================================
axins3 = inset_axes(
    ax3,
    width="34%",
    height="34%",
    loc='center right',
    borderpad=1
)

axins3.plot(Jof4, Id4, color='purple', lw=1.8)
axins3.plot(Jof9, Id9, color='black', lw=1.8)
axins3.plot(Jof5, Id5, color='brown', lw=1.8)

axins3.set_xscale("log")

axins3.set_xlim(9e-3, 2)
axins3.set_ylim(0.0042, 0.00431)

axins3.tick_params(direction='in', labelsize=5)

for spine in axins3.spines.values():
    spine.set_linewidth(0.6)

mark_inset(
    ax3,
    axins3,
    loc1=2,
    loc2=4,
    fc="none",
    ec="0.4",
    lw=0.6
)

# ==========================================================
# (d) EFFICIENCY
# ==========================================================
datasets_d = [eff0, eff2, eff3, eff9, eff4, eff5, eff14]

for (x, c, lab), y in zip(curves, datasets_d):
    ls = '--' if '10^{2}' in lab else '-'
    ax4.plot(x, y, color=c, lw=LINE_W, ls=ls)

ax4.set_xscale('log')

ax4.set_ylabel(
    r'$\eta_2$',
    fontsize=LABEL_FS
)

ax4.text(
    0.90, 0.89, '(d)',
    transform=ax4.transAxes,
    fontsize=PANEL_FS,
    fontweight='bold'
)

# ==========================================================
# INSET (d)
# ==========================================================
axins4 = inset_axes(
    ax4,
    width="34%",
    height="34%",
    loc='center right',
    borderpad=1
)

axins4.plot(Jof4, eff4, color='purple', lw=1.8)
axins4.plot(Jof9, eff9, color='black', lw=1.8)
axins4.plot(Jof5, eff5, color='brown', lw=1.8)

axins4.set_xscale("log")

axins4.set_xlim(1e-2, 4)
axins4.set_ylim(0.5, 0.541)

axins4.tick_params(direction='in', labelsize=5)

for spine in axins4.spines.values():
    spine.set_linewidth(0.6)

mark_inset(
    ax4,
    axins4,
    loc1=2,
    loc2=4,
    fc="none",
    ec="0.4",
    lw=0.6
)

# ==========================================================
# SHARED FORMATTING
# ==========================================================
for ax in axs.flat:

    ax.tick_params(
        direction='in',
        which='both',
        labelsize=TICK_FS
    )

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

# ==========================================================
# SHARED X LABELS
# ==========================================================
ax3.set_xlabel(
    r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_L)$',
    fontsize=LABEL_FS
)

ax4.set_xlabel(
    r'$J_{0}/(\beta_{\mathrm{ph}}\kappa_L)$',
    fontsize=LABEL_FS
)

# ==========================================================
# GLOBAL LEGEND
# ==========================================================
handles, labels = ax1.get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc='upper center',
    ncol=4,
    fontsize=LEG_FS,
    frameon=True
)

# ==========================================================
# LAYOUT
# ==========================================================
plt.tight_layout(
    pad=0.6,
    rect=[0, 0, 1, 0.92]
)

# ==========================================================
# SAVE
# ==========================================================
plt.savefig(
    'fig_PRB_2x2new.pdf',
    bbox_inches='tight'
)

plt.show()