import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


data0 = np.load("phonong=0b100.npz")
Jof0,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0,Act0,Nqm0,Nltotal0,Work0,eff0,effph0,coheveig0 = data0["Jof"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["concv"], data0["Nls"],data0["Acts"], data0["Nlqm"],data0["Nltotal"],data0["Work"],data0["eff"],data0["effph"],data0["coheveig"]

data1 = np.load("phonong=10^{-4}b100.npz")
Jof1,Id1,Ile1,Ire1,Iphs1,cohes1,concv1,Nls1,Act1,Nqm1,Nltotal1,Work1,eff1,effph1,coheveig1 = data1["Jof"], data1["Id"], data1["Ile"], data1["Ire"], data1["Iphs"], data1["cohes"], data1["concv"], data1["Nls"],data1["Acts"], data1["Nlqm"],data1["Nltotal"],data1["Work"],data1["eff"],data1["effph"],data1["coheveig"]   

data2 = np.load("phonong=5_10^{-4}b100.npz")
Jof2,Id2,Ile2,Ire2,Iphs2,cohes2,concv2,Nls2,Act2,Nqm2,Nltotal2,Work2,eff2,effph2,coheveig2 = data2["Jof"], data2["Id"], data2["Ile"], data2["Ire"], data2["Iphs"], data2["cohes"], data2["concv"], data2["Nls"],data2["Acts"], data2["Nlqm"],data2["Nltotal"],data2["Work"],data2["eff"],data2["effph"],data2["coheveig"]
data3 = np.load("phonong=10^{-3}b100.npz")
Jof3,Id3,Ile3,Ire3,Iphs3,cohes3,concv3,Nls3,Act3,Nqm3,Nltotal3,Work3,eff3,effph3,coheveig3 = data3["Jof"], data3["Id"], data3["Ile"], data3["Ire"], data3["Iphs"], data3["cohes"], data3["concv"], data3["Nls"],data3["Acts"], data3["Nlqm"],data3["Nltotal"],data3["Work"],data3["eff"],data3["effph"],data3["coheveig"]

data4 = np.load("phonong=5_10^{-3}b100.npz")
Jof4,Id4,Ile4,Ire4,Iphs4,cohes4,concv4,Nls4,Act4,Nqm4,Nltotal4,Work4,eff4,effph4,coheveig4 = data4["Jof"], data4["Id"], data4["Ile"], data4["Ire"], data4["Iphs"], data4["cohes"], data4["concv"], data4["Nls"],data4["Acts"], data4["Nlqm"],data4["Nltotal"],data4["Work"],data4["eff"],data4["effph"],data4["coheveig"]

data5 = np.load("phonong=10^{-2}b100.npz")
Jof5,Id5,Ile5,Ire5,Iphs5,cohes5,concv5,Nls5,Act5,Nqm5,Nltotal5,Work5,eff5,effph5,coheveig5 = data5["Jof"], data5["Id"], data5["Ile"], data5["Ire"], data5["Iphs"], data5["cohes"], data5["concv"], data5["Nls"],data5["Acts"], data5["Nlqm"],data5["Nltotal"],data5["Work"],data5["eff"],data5["effph"],data5["coheveig"]

data6 = np.load("phonong=5_10^{-2}b100.npz")
Jof6,Id6,Ile6,Ire6,Iphs6,cohes6,concv6,Nls6,Act6,Nqm6,Nltotal6,Work6,eff6,effph6,coheveig6 = data6["Jof"], data6["Id"], data6["Ile"], data6["Ire"], data6["Iphs"], data6["cohes"], data6["concv"], data6["Nls"],data6["Acts"], data6["Nlqm"],data6["Nltotal"],data6["Work"],data6["eff"],data6["effph"],data6["coheveig"]


data7 = np.load("phonong=10^{-1}b100.npz")
Jof7,Id7,Ile7,Ire7,Iphs7,cohes7,concv7,Nls7,Act7,Nqm7,Nltotal7,Work7,eff7,effph7,coheveig7 = data7["Jof"], data7["Id"], data7["Ile"], data7["Ire"], data7["Iphs"], data7["cohes"], data7["concv"], data7["Nls"],data7["Acts"], data7["Nlqm"],data7["Nltotal"],data7["Work"],data7["eff"],data7["effph"],data7["coheveig"]


data8 = np.load("phonong=7_10^{-3}b100.npz")
Jof8,Id8,Ile8,Ire8,Iphs8,cohes8,concv8,Nls8,Act8,Nqm8,Nltotal8,Work8,eff8,effph8,coheveig8 = data8["Jof"], data8["Id"], data8["Ile"], data8["Ire"], data8["Iphs"], data8["cohes"], data8["concv"], data8["Nls"],data8["Acts"], data8["Nlqm"],data8["Nltotal"],data8["Work"],data8["eff"],data8["effph"],data8["coheveig"]


data9 = np.load("phonong=3_10^{-3}b100.npz")
Jof9,Id9,Ile9,Ire9,Iphs9,cohes9,concv9,Nls9,Act9,Nqm9,Nltotal9,Work9,eff9,effph9,coheveig9 = data9["Jof"], data9["Id"], data9["Ile"], data9["Ire"], data9["Iphs"], data9["cohes"], data9["concv"], data9["Nls"],data9["Acts"], data9["Nlqm"],data9["Nltotal"],data9["Work"],data9["eff"],data9["effph"],data9["coheveig"]

data10 = np.load("phonong=10b100.npz")
Jof10,Id10,Ile10,Ire10,Iphs10,cohes10,concv10,Nls10,Work10,eff10,effph10,coheveig10 = data10["Jof"], data10["Id"], data10["Ile"], data10["Ire"], data10["Iphs"], data10["cohes"], data10["concv"], data10["Nls"],data10["Work"],data10["eff"],data10["effph"],data10["coheveig"]    

data11 = np.load("phonong=100b100.npz")
Jof11,Id11,Ile11,Ire11,Iphs11,cohes11,concv11,Nls11,Work11,eff11,effph11,coheveig11 = data11["Jof"], data11["Id"], data11["Ile"], data11["Ire"], data11["Iphs"], data11["cohes"], data11["concv"], data11["Nls"],data11["Work"],data11["eff"],data11["effph"],data11["coheveig"]

data12 = np.load("phonong=1000b100.npz")
Jof12,Id12,Ile12,Ire12,Iphs12,cohes12,concv12,Nls12,Work12,eff12,effph12,coheveig12 = data12["Jof"], data12["Id"], data12["Ile"], data12["Ire"], data12["Iphs"], data12["cohes"], data12["concv"], data12["Nls"],data12["Work"],data12["eff"],data12["effph"],data12["coheveig"]

data13 = np.load("phonong=10000b100.npz")
Jof13,Id13,Ile13,Ire13,Iphs13,cohes13,concv13,Nls13,Work13,eff13,effph13,coheveig13 = data13["Jof"], data13["Id"], data13["Ile"], data13["Ire"], data13["Iphs"], data13["cohes"], data13["concv"], data13["Nls"],data13["Work"],data13["eff"],data13["effph"],data13["coheveig"]

data14 = np.load("phonong=1b100.npz")  
Jof14,Id14,Ile14,Ire14,Iphs14,cohes14,concv14,Nls14,Work14,eff14,effph14,coheveig14 = data14["Jof"], data14["Id"], data14["Ile"], data14["Ire"], data14["Iphs"], data14["cohes"], data14["concv"], data14["Nls"],data14["Work"],data14["eff"],data14["effph"],data14["coheveig"]

data15 = np.load("phonong=5b100.npz")
Jof15,Id15,Ile15,Ire15,Iphs15,cohes15,concv15,Nls15,Work15,eff15,effph15,coheveig15 = data15["Jof"], data15["Id"], data15["Ile"], data15["Ire"], data15["Iphs"], data15["cohes"], data15["concv"], data15["Nls"],data15["Work"],data15["eff"],data15["effph"],data15["coheveig"]

data16 = np.load("phonong=9b100.npz")
Jof16,Id16,Ile16,Ire16,Iphs16,cohes16,concv16,Nls16,Work16,eff16,effph16,coheveig16 = data16["Jof"], data16["Id"], data16["Ile"], data16["Ire"], data16["Iphs"], data16["cohes"], data16["concv"], data16["Nls"],data16["Work"],data16["eff"],data16["effph"],data16["coheveig"]
######################
#####seculardata######
######################


N00 = len(Nls0)

for i in range(N00):
    Nls0[i] = Nls0[i]*100
    Nls1[i] = Nls1[i]*100
    Nls2[i] = Nls2[i]*100
    Nls3[i] = Nls3[i]*100
    Nls4[i] = Nls4[i]*100
    Nls5[i] = Nls5[i]*100
    Nls6[i] = Nls6[i]*100
    Nls7[i] = Nls7[i]*100
    Nls8[i] = Nls8[i]*100
    Nls9[i] = Nls9[i]*100
    Nls10[i] = Nls10[i]*100
    Nls11[i] = Nls11[i]*100
    Nls12[i] = Nls12[i]*100
    Nls13[i] = Nls13[i]*100
    Nls14[i] = Nls14[i]*100
    Nls15[i] = Nls15[i]*100
    Nls16[i] = Nls16[i]*100



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

n = len(Jof0)
Td = 2
T = 100
effcarnot = []

for i in range(n):
    eta = 1-(Td/T)
    effcarnot.append(eta)



plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 


plt.plot(Jof0,eff0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,eff1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,eff2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,eff3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof9,eff9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.plot(Jof4,eff4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,eff5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
#plt.plot(Jof6,eff6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,eff7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
#plt.plot(Jof8,eff8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
#plt.plot(Jof0,effcarnot,color='blue',lw=3,linestyle = '--', label = r'$\eta_{carnot}$')
#plt.xlabel(r'$J_{0}/\beta_{Ph}\kappa_{L}$',fontsize = 28)
#plt.ylabel(r'$\eta_{LR}$',fontsize = 28)
plt.xticks(fontsize=30)  # X-axis tick labels
plt.yticks(fontsize=30)
#plt.legend(fontsize=22,loc = "lower right")
#plt.figtext(0.14, 0.955, '(a)',fontsize=30, fontweight='bold', va='top', ha='left')

plt.xscale("log")
plt.show()

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
ax1.plot(Jof0, eff0, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=0$')
#ax1.plot(Jof1, eff1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
ax1.plot(Jof2, eff2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\times 10^{-2}$')
ax1.plot(Jof3, eff3, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{-1}$')
ax1.plot(Jof9, eff9, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=3\times10^{-1}$')
ax1.plot(Jof4, eff4, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=5\times10^{-1}$')
ax1.plot(Jof5, eff5, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{0}$')
#ax1.plot(Jof7, eff7, color='gray',   lw=LINE_W, label=r'$g/\kappa_{L}=10^{1}$')
#ax1.plot(Jof10, eff10, color='red',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{3}$')
#ax1.plot(Jof11, eff11, color='blue',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{4}$')
#ax1.plot(Jof12, eff12, color='green',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{5}$')
#ax1.plot(Jof13, eff13, color='brown',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{6}$')


###########zoomplot###############
'''x1 = Jof4[10:40]
y1 = eff4[10:40]

zm = ax1.inset_axes([0.4,0.5,0.5,0.5]) #x0,y0,width height
zm.plot(x1,y1)
ax1.indicate_inset_zoom(zm,edgecolor = 'blue')'''

axins = inset_axes(
    ax1,
    width="50%",
    height="50%",
    bbox_to_anchor=(0.35, 0.05, 0.6, 0.6),
    bbox_transform=ax1.transAxes,
    
    borderpad=0
)

# Plot SAME curves (no legend!)
axins.plot(Jof4, eff4, color='purple', lw=2)
axins.plot(Jof9, eff9, color='black',  lw=2)

# Zoom region (example values — adjust!)
axins.set_xlim(7e-4, 1e-0)
axins.set_ylim(0.5, 0.535)

# Inset formatting
axins.tick_params(direction='in', labelsize=6)
for spine in axins.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)


ax1.set_xscale("log")
ax1.set_ylabel(r'$\eta_{2}$', fontsize=LABEL_FS)

ax1.text(0.92, 0.90, '(a)', transform=ax1.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
ax2.plot(Jof0, cohes0, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, cohes1, color='orange', lw=LINE_W)
ax2.plot(Jof2, cohes2, color='green',  lw=LINE_W)
ax2.plot(Jof3, cohes3, color='red',    lw=LINE_W)
ax2.plot(Jof9, cohes9, color='black',  lw=LINE_W)
ax2.plot(Jof4, cohes4, color='purple', lw=LINE_W)
ax2.plot(Jof5, cohes5, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, cohes7, color='gray',   lw=LINE_W)
#ax2.plot(Jof10, cohes10, color='red',   lw=LINE_W, ls = '--')
#ax2.plot(Jof11, cohes11, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, cohes12, color='green',   lw=LINE_W, ls = '--')
#ax2.plot(Jof13, cohes13, color='brown',   lw=LINE_W, ls = '--')



ax2.set_xscale("log")
ax2.set_ylabel(r'$\mathcal{C}_{l_{1}}$', fontsize=LABEL_FS)

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
    ncol=3,
    fontsize=LEG_FS,
    frameon=True
)

# ===================== Layout & save =====================
plt.tight_layout(pad=0.4, rect=[0, 0, 1, 0.88])
plt.savefig("figeffcohe1.pdf")
plt.close()




'''
plt.plot(Jof0,Work0, color='blue',lw=3, label = r'$g/\kappa_{L}=0$')
plt.plot(Jof1,Work1, color='orange',lw=3, label = r'$g/\kappa_{L}=10^{-2}$')
plt.plot(Jof2,Work2, color='green',lw=3, label = r'$g/\kappa_{L}=5\cdot10^{-2}$')
plt.plot(Jof3,Work3, color='red',lw=3, label = r'$g/\kappa_{L}=10^{-1}$')
plt.plot(Jof9,Work9, color='black',lw=3,linestyle = '--', label = r'$g/\kappa_{L}=3\cdot 10^{-1}$')
plt.plot(Jof4,Work4, color='purple',lw=3, label = r'$g/\kappa_{L}=5\cdot10^{-1}$')
plt.plot(Jof5,Work5, color='brown',lw=3, label = r'$g/\kappa_{L}=10^{0}$')
#plt.plot(Jof6,Work6, color='pink',lw=3, label = r'$g/\kappa_{L}=5\cdot10^{0}$')
plt.plot(Jof7,Work7, color='gray',lw=3, label = r'$g/\kappa_{L}=10^{1}$')
#plt.plot(Jof8,Work8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')
#plt.plot(Jof9,Work9, color='black',lw=3,linestyle = '--', label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.xlabel(r'$J_{0}/\beta_{Ph}\kappa_{L}$',fontsize = 28)
plt.ylabel(r'$\dot{W}_{LR}$',fontsize = 28)
plt.xticks(fontsize=30)  # X-axis tick labels
plt.yticks(fontsize=30)
plt.legend(fontsize=22,loc = "upper right")
plt.figtext(0.14, 0.955, '(a)',fontsize=30, fontweight='bold', va='top', ha='left')

plt.xscale("log")
plt.show()



plt.plot(Jof0,cohes0, color='blue',lw=3, label = r'$\frac{g}{\kappa_{L}}=0$')
plt.plot(Jof1,cohes1, color='orange',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-2}$')
plt.plot(Jof2,cohes2, color='green',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-2}$')
plt.plot(Jof3,cohes3, color='red',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{-1}$')
plt.plot(Jof9,cohes9, color='black',linestyle = '--',lw=3, label = r'$\frac{g}{\kappa_{L}}=3\cdot 10^{-1}$')
plt.plot(Jof4,cohes4, color='purple',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{-1}$')
plt.plot(Jof5,cohes5, color='brown',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{0}$')
#plt.plot(Jof6,cohes6, color='pink',lw=3, label = r'$\frac{g}{\kappa_{L}}=5\cdot10^{0}$')
plt.plot(Jof7,cohes7, color='gray',lw=3, label = r'$\frac{g}{\kappa_{L}}=10^{1}$')
#plt.plot(Jof8,cohes8, color='black',lw=3, label = r'$\frac{g}{\kappa_{L}}=7\cdot 10^{-1}$')

plt.xlabel(r'$J_{0}/\beta_{Ph}\kappa_{L}$',fontsize = 28)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$',fontsize = 28)
plt.xticks(fontsize=30)  # X-axis tick labels
plt.yticks(fontsize=30)
plt.figtext(0.14, 0.955, '(b)',fontsize=30, fontweight='bold', va='top', ha='left')

#plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()'''


datos1 = np.load("phonong=0.1b100sec.npz")
Jof1s,Id1s,Ile1s,Ire1s,Iphs1s,cohes1s,concv1s,Nls1s,Work1s,eff1s,effph1s = datos1["Jof"], datos1["Id"], datos1["Ile"], datos1["Ire"], datos1["Iphs"], datos1["cohes"], datos1["concv"], datos1["Nls"],datos1["Work"],datos1["eff"],datos1["effph"]

datos2 = np.load("phonong=1b100sec.npz")
Jof2s,Id2s,Ile2s,Ire2s,Iphs2s,cohes2s,concv2s,Nls2s,Work2s,eff2s,effph2s = datos2["Jof"], datos2["Id"], datos2["Ile"], datos2["Ire"], datos2["Iphs"], datos2["cohes"], datos2["concv"], datos2["Nls"],datos2["Work"],datos2["eff"],datos2["effph"]

datos3 = np.load("phonong=10b100sec.npz")
Jof3s,Id3s,Ile3s,Ire3s,Iphs3s,cohes3s,concv3s,Nls3s,Work3s,eff3s,effph3s = datos3["Jof"], datos3["Id"], datos3["Ile"], datos3["Ire"], datos3["Iphs"], datos3["cohes"], datos3["concv"], datos3["Nls"],datos3["Work"],datos3["eff"],datos3["effph"]

datos4 = np.load("phonong=100b100sec.npz")
Jof4s,Id4s,Ile4s,Ire4s,Iphs4s,cohes4s,concv4s,Nls4s,Work4s,eff4s,effph4s = datos4["Jof"], datos4["Id"], datos4["Ile"], datos4["Ire"], datos4["Iphs"], datos4["cohes"], datos4["concv"], datos4["Nls"],datos4["Work"],datos4["eff"],datos4["effph"]

datos5 = np.load("phonong=1000b100sec.npz")
Jof5s,Id5s,Ile5s,Ire5s,Iphs5s,cohes5s,concv5s,Nls5s,Work5s,eff5s,effph5s = datos5["Jof"], datos5["Id"], datos5["Ile"], datos5["Ire"], datos5["Iphs"], datos5["cohes"], datos5["concv"], datos5["Nls"],datos5["Work"],datos5["eff"],datos5["effph"]



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
#ax1.plot(Jof1s, eff1s, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=10$')
#ax1.plot(Jof1, eff1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
#ax1.plot(Jof2, eff2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\cdot10^{-2}$')
ax1.plot(Jof2s, eff2s, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{2}$')
ax1.plot(Jof3s, eff3s, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{3}$', ls = '--')
ax1.plot(Jof4s, eff4s, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=10^{4}$')
ax1.plot(Jof5s, eff5s, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{5}$')


ax1.set_xscale("log")
ax1.set_ylabel(r'$\eta_{LR}$', fontsize=LABEL_FS)

ax1.text(0.92, 0.90, '(a)', transform=ax1.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
#ax2.plot(Jof1s, cohes1s, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, cohes1, color='orange', lw=LINE_W)
#ax2.plot(Jof2, cohes2, color='green',  lw=LINE_W)
ax2.plot(Jof2s, cohes2s, color='red',    lw=LINE_W)
ax2.plot(Jof3s, cohes3s, color='black',  lw=LINE_W, ls='--')
ax2.plot(Jof4s, cohes4s, color='purple', lw=LINE_W)
ax2.plot(Jof5s, cohes5s, color='brown',  lw=LINE_W)


ax2.set_xscale("log")
ax2.set_ylabel(r'$\mathcal{C}_{l_{1}}$', fontsize=LABEL_FS)

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
plt.savefig("figeffcohesec.pdf")
plt.show()
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
ax1.plot(Jof0, eff0, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=0$')
#ax1.plot(Jof1, eff1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
ax1.plot(Jof2, eff2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\times 10^{-2}$')
ax1.plot(Jof3, eff3, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{-1}$')
ax1.plot(Jof9, eff9, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=3\times10^{-1}$')
ax1.plot(Jof4, eff4, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=5\times10^{-1}$')
ax1.plot(Jof5, eff5, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{0}$')
#ax1.plot(Jof7, eff7, color='gray',   lw=LINE_W, label=r'$g/\kappa_{L}=10^{1}$')
#ax1.plot(Jof10, eff10, color='red',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{3}$')
#ax1.plot(Jof11, eff11, color='blue',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{4}$')
#ax1.plot(Jof12, eff12, color='green',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{5}$')
#ax1.plot(Jof13, eff13, color='brown',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{6}$')
ax1.plot(Jof14, eff14, color='black',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{2}$')
#ax1.plot(Jof15, eff15, color='gray',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=5\times 10^{2}$')
#ax1.plot(Jof16, eff16, color='pink',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=9\times 10^{2}$')

axins = inset_axes(
    ax1,
    width="50%",
    height="50%",
    bbox_to_anchor=(0.35, 0.05, 0.6, 0.6),
    bbox_transform=ax1.transAxes,
    
    borderpad=0
)

# Plot SAME curves (no legend!)
axins.plot(Jof4, eff4, color='purple', lw=2)
axins.plot(Jof5,eff5,color='brown', lw = 2)
axins.plot(Jof9, eff9, color='black',  lw=2)

# Zoom region (example values — adjust!)
axins.set_xlim(7e-4, 2.1e-0)
axins.set_ylim(0.5, 0.542)

# Inset formatting
axins.tick_params(direction='in', labelsize=6)
for spine in axins.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)


ax1.set_xscale("log")
ax1.set_ylabel(r'$\eta_{2}$', fontsize=LABEL_FS)

ax1.text(0.92, 0.90, '(a)', transform=ax1.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
ax2.plot(Jof0, coheveig0, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, coheveig1, color='orange', lw=LINE_W)
ax2.plot(Jof2, coheveig2, color='green',  lw=LINE_W)
ax2.plot(Jof3, coheveig3, color='red',    lw=LINE_W)
ax2.plot(Jof9, coheveig9, color='black',  lw=LINE_W)
ax2.plot(Jof4, coheveig4, color='purple', lw=LINE_W)
ax2.plot(Jof5, coheveig5, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, coheveig7, color='gray',   lw=LINE_W)
#ax2.plot(Jof10, coheveig10, color='red',   lw=LINE_W, ls = '--')
#ax2.plot(Jof11, coheveig11, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, coheveig12, color='green',   lw=LINE_W, ls = '--')
ax2.plot(Jof13, coheveig13, color='brown',   lw=LINE_W, ls = '--')
ax2.plot(Jof14, coheveig14, color='black',   lw=LINE_W, ls = '--')
#ax2.plot(Jof15, coheveig15, color='gray',   lw=LINE_W, ls = '--')
#ax2.plot(Jof16, coheveig16, color='pink',   lw=LINE_W, ls = '--')



ax2.set_xscale("log")
ax2.set_ylabel(r'$\mathcal{C}_{l_{1}}(\mathrm{eig})$', fontsize=LABEL_FS)

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
plt.savefig("figeffcoheeig.pdf")
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
    width="32%",
    height="32%",
    #bbox_to_anchor=(0.005, 0.005, 0.1, 0.2),  #[x, y, width, height] 
    #bbox_transform=ax10.transAxes,
    loc = 'lower right',
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
ax20.plot(Jof0, Work0, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, coheveig1, color='orange', lw=LINE_W)
ax20.plot(Jof2, Work2, color='green',  lw=LINE_W)
ax20.plot(Jof3, Work3, color='red',    lw=LINE_W)
ax20.plot(Jof9, Work9, color='black',  lw=LINE_W)
ax20.plot(Jof4, Work4, color='purple', lw=LINE_W)
ax20.plot(Jof5, Work5, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, coheveig7, color='gray',   lw=LINE_W)
#ax20.plot(Jof10, Work10, color='red',   lw=LINE_W, ls = '--')
#ax20.plot(Jof11, Work11, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, coheveig12, color='green',   lw=LINE_W, ls = '--')
#ax20.plot(Jof13, Work13, color='brown',   lw=LINE_W, ls = '--')
ax20.plot(Jof14, Work14, color='black',   lw=LINE_W, ls = '--')
#ax2.plot(Jof15, coheveig15, color='gray',   lw=LINE_W, ls = '--')
#ax2.plot(Jof16, coheveig16, color='pink',   lw=LINE_W, ls = '--')



ax20.set_xscale("log")
ax20.set_ylabel(r'$\dot{W}_{2}$', fontsize=LABEL_FS)

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
axins20.plot(Jof4, Work4, color='purple', lw=2)
axins20.plot(Jof9, Work9, color='black',  lw=2)
axins20.plot(Jof5, Work5, color='brown',  lw=2)
axins20.set_xscale("log")

# Zoom region (example values — adjust!)
axins20.set_xlim(1e-2, 4)
axins20.set_ylim(-0.048, -0.039)

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
plt.savefig("figinfocoheeig.pdf")
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
ax1.plot(Jof0, Nls0, color='blue',   lw=LINE_W, label=r'$g/\kappa_{L}=0$')
#ax1.plot(Jof1, Nls1, color='orange', lw=LINE_W, label=r'$g/\kappa_{L}=10^{-2}$')
ax1.plot(Jof2, Nls2, color='green',  lw=LINE_W, label=r'$g/\kappa_{L}=5\times 10^{-2}$')
ax1.plot(Jof3, Nls3, color='red',    lw=LINE_W, label=r'$g/\kappa_{L}=10^{-1}$')
ax1.plot(Jof9, Nls9, color='black',  lw=LINE_W, label=r'$g/\kappa_{L}=3\times10^{-1}$')
ax1.plot(Jof4, Nls4, color='purple', lw=LINE_W, label=r'$g/\kappa_{L}=5\times10^{-1}$')
ax1.plot(Jof5, Nls5, color='brown',  lw=LINE_W, label=r'$g/\kappa_{L}=10^{0}$')
#ax1.plot(Jof7, Nls7, color='gray',   lw=LINE_W, label=r'$g/\kappa_{L}=10^{1}$')
#ax1.plot(Jof10, Nls10, color='red',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{3}$')
#ax1.plot(Jof11, Nls11, color='blue',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{4}$')
#ax1.plot(Jof12, Nls12, color='green',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{5}$')
#ax1.plot(Jof13, Nls13, color='brown',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{6}$')
ax1.plot(Jof14, Nls14, color='black',   lw=LINE_W, ls = '--',label=r'$g/\kappa_{L}=10^{2}$')
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
axins.set_xlim(7e-4, 2.1e-0)
axins.set_ylim(-0.05, -0.0384)

# Inset formatting
axins.tick_params(direction='in', labelsize=6)
for spine in axins.spines.values():
    spine.set_linewidth(0.6)

# Optional: draw connectors
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.4", lw=0.6)


ax1.set_xscale("log")
ax1.set_ylabel(r'$\dot{N}_{B_L}/\kappa_{L}$', fontsize=LABEL_FS)

ax1.text(0.92, 0.90, '(a)', transform=ax1.transAxes,
         fontsize=PANEL_FS, fontweight='bold')

# ===================== (b) Coherence =====================
ax2.plot(Jof0, coheveig0, color='blue',   lw=LINE_W)
#ax2.plot(Jof1, coheveig1, color='orange', lw=LINE_W)
ax2.plot(Jof2, coheveig2, color='green',  lw=LINE_W)
ax2.plot(Jof3, coheveig3, color='red',    lw=LINE_W)
ax2.plot(Jof9, coheveig9, color='black',  lw=LINE_W)
ax2.plot(Jof4, coheveig4, color='purple', lw=LINE_W)
ax2.plot(Jof5, coheveig5, color='brown',  lw=LINE_W)
#ax2.plot(Jof7, coheveig7, color='gray',   lw=LINE_W)
#ax2.plot(Jof10, coheveig10, color='red',   lw=LINE_W, ls = '--')
#ax2.plot(Jof11, coheveig11, color='blue',   lw=LINE_W, ls = '--')
#ax2.plot(Jof12, coheveig12, color='green',   lw=LINE_W, ls = '--')
ax2.plot(Jof13, coheveig13, color='brown',   lw=LINE_W, ls = '--')
ax2.plot(Jof14, coheveig14, color='black',   lw=LINE_W, ls = '--')
#ax2.plot(Jof15, coheveig15, color='gray',   lw=LINE_W, ls = '--')
#ax2.plot(Jof16, coheveig16, color='pink',   lw=LINE_W, ls = '--')



ax2.set_xscale("log")
ax2.set_ylabel(r'$\mathcal{C}_{l_{1}}(\mathrm{eig})$', fontsize=LABEL_FS)

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
plt.savefig("figcurrenteignew.pdf")
plt.close()
