import numpy as np
import matplotlib.pyplot as plt

data0 = np.load("phononJ0=10^{-1}gpaper.npz")
gof10,Qlr0,Qph0,Id0,Ile0,Ire0,Iphs0,cohes0,Ilrnew0,Eds0,Els0,Ers0,Work0,concv0,Nls0,coheveig0 = data0["gof1"], data0["Qlr"], data0["Qphs"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["Ilrnew"], data0["Eds"], data0["Els"], data0["Ers"], data0["Work"], data0["concv"], data0["Nls"], data0["coheveig"]

data1 = np.load("phononJ0=10^{-2}gpaper.npz")
gof11,Qlr1,Qph1,Id1,Ile1,Ire1,Iphs1,cohes1,Ilrnew1,Eds1,Els1,Ers1,Work1,concv1,Nls1,coheveig1 = data1["gof1"], data1["Qlr"], data1["Qphs"], data1["Id"], data1["Ile"], data1["Ire"], data1["Iphs"], data1["cohes"], data1["Ilrnew"], data1["Eds"], data1["Els"], data1["Ers"], data1["Work"], data1["concv"], data1["Nls"], data1["coheveig"]

data2 = np.load("phononJ0=0gpaper.npz")
gof12,Qlr2,Qph2,Id2,Ile2,Ire2,Iphs2,cohes2,Ilrnew2,Eds2,Els2,Ers2,Work2,concv2,Nls2,coheveig2 = data2["gof1"], data2["Qlr"], data2["Qphs"], data2["Id"], data2["Ile"], data2["Ire"], data2["Iphs"], data2["cohes"], data2["Ilrnew"], data2["Eds"], data2["Els"], data2["Ers"], data2["Work"], data2["concv"], data2["Nls"], data2["coheveig"]
T = 100

Ti0 = []
Ti1 = []
Ti2 = []

f0 =[]
f1 = []
f2 = []
E22s = []
E21s = []
E20s = []

Wds = []

for i in range(len(gof10)):
    E22 = Els2[i] +Ers2[i]
    E21 = Els1[i] +Ers1[i]
    E20 = Els0[i] +Ers0[i]
    Ti0.append(T*Ile0[i])
    Ti1.append(T*Ile1[i])
    Ti2.append(T*Ile2[i])
    f0.append(E20+ T*Ile0[i])
    f1.append(E21+ T*Ile1[i])
    f2.append(E22+ T*Ile2[i])
    E22s.append(E22)
    E21s.append(E21)
    E20s.append(E20)
    Wds.append(0)
# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns


#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax10.plot(gof10,Iphs0, color='red',lw=3, label = r'$\dot{I}_{ph}$')
ax10.plot(gof10,Ile0, color='black',lw=3, label = r'$\dot{I}_{Le}$')
ax10.plot(gof10,Ire0, color='blue',lw=3, label = r'$\dot{I}_{Re}$')
ax10.set_ylabel(r'$\dot{I}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "center left")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


ax20.plot(gof10,cohes0,label = r'$\mathcal{C}_{l_{1}}$', color = 'b',lw = 3)
ax20.plot(gof10,concv0, label = r'$\mathcal{C}_{on}$', color = 'r',lw=3)  
ax20.set_xlabel(r'$g_{0}/\kappa_{L}$',fontsize = 20)   
ax20.set_xscale("log")
ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia y entrelazamiento",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')

plt.tight_layout()  # Avoids overlapping labels
plt.show()


###mejorar formato############
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# Create subplots (1 row, 2 columns)
fig, (ax11, ax21) = plt.subplots(2, 1,sharex=True, figsize=(3.39, 4.5))  # 1 row, 2 columns

ax11.plot(gof12,Nls2,color='black',lw=LINE_W, label = r'$\dot{N}_{B_L}/\kappa_{L}$')
#ax11.plot(gof10,Work0, color='red',linestyle = '--',lw=LINE_W)
#(0.7,0.48)
#ax11.set_ylabel(r'$I/\kappa_{L}$',fontsize = 20)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
#ax10.xticks(fontsize=17)  
#ax10.yticks(fontsize=17)
#ax11.legend(bbox_to_anchor=(0.20, 0.98),fontsize=15,loc = "upper left", ncol = 2)
#ax11.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax11.tick_params(labelbottom=False,direction='in', which='both',labelsize = TICK_FS)
ax11.text(0.90, 0.94, '(a)', transform=ax11.transAxes, fontsize=PANEL_FS, fontweight='bold')
ax11.set_xscale("log")
ax11.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left',
    bbox_to_anchor=(0.7, 0.25)
)

ax21.plot(gof12,cohes2,label = r'$\mathcal{C}_{l_{1}}(\mathrm{eig})$', color = 'black',lw = LINE_W)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
ax21.plot(gof12,concv2, label = r'$\mathcal{C}_{on}$', color = 'r',lw=LINE_W) 
#plt.plot(eVs,Qdf,label = r'$J_{d}$',color = "gray",lw=2)
#ax2.xticks(fontsize=17)  # X-axis tick labels
#ax2.yticks(fontsize=17)  # Y-axis tick labels
#plt.xscale("log")
ax21.set_xlabel(r'$g/\kappa_{L}$',fontsize = LABEL_FS)
#ax21.axvspan(0, 2.465, facecolor='b', alpha=0.5)
#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  

ax21.tick_params(labelsize=TICK_FS)  # font size of tick labels 
ax21.text(0.90, 0.94, '(b)', transform=ax21.transAxes, fontsize=PANEL_FS, fontweight='bold')

ax21.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left',
    bbox_to_anchor=(0.7, 0.25)
)

ax21.set_xscale("log")
#fig.supylabel("Cantidades termodinámicas", fontsize=22)
#plt.subplots_adjust(left=0.05) 
plt.tight_layout(pad=0.4)
plt.savefig("figcurrentg.pdf")
plt.show()
plt.close()


# Create subplots (1 row, 2 columns)
fig, (ax11, ax21) = plt.subplots(2, 1,sharex=True, figsize=(3.39, 4.5))  # 1 row, 2 columns

ax11.plot(gof10,E22s,color='blue',lw=LINE_W, label = r'$\dot{E}_{2}$')
ax11.plot(gof10,Work2,color='red',lw=LINE_W, label = r'$\dot{W}_{2}$')
ax11.plot(gof10,Ti2,color='black',linestyle = '--',lw=LINE_W, label = r'$T\dot{I}_{2}$')
ax11.plot(gof10,f2,color='black',lw=LINE_W, label = r'$\dot{\mathcal{F}}_{2}$')
ax11.plot(gof10,Qlr2,color='purple',lw=LINE_W, label = r'$\dot{Q}_{2}$')
#(0.7,0.48)
#ax11.set_ylabel(r'$I/\kappa_{L}$',fontsize = 20)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
#ax10.xticks(fontsize=17)  
#ax10.yticks(fontsize=17)
#ax11.legend(bbox_to_anchor=(0.20, 0.98),fontsize=15,loc = "upper left", ncol = 2)
#ax11.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax11.tick_params(labelbottom=False,direction='in', which='both',labelsize = TICK_FS)
ax11.text(0.10, 0.89, '(a)', transform=ax11.transAxes, fontsize=PANEL_FS, fontweight='bold')
ax11.set_xscale("log")
ax11.legend(
    fontsize=7,
    frameon=True,
    ncol=2,
    loc='upper right',
    bbox_to_anchor=(0.43, 0.30)
)

ax21.plot(gof10,Eds2,label = r'$\dot{E}_{1}$', color = 'blue',lw = LINE_W)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
ax21.plot(gof10,Wds,label = r'$\dot{W}_{1}$', color = 'red',lw=LINE_W)
#ax21.plot(gof10,concv2,label = r'$\mathcal{C}_{on}$', color = 'r',lw=LINE_W) 
#plt.plot(eVs,Qdf,label = r'$J_{d}$',color = "gray",lw=2)
#ax2.xticks(fontsize=17)  # X-axis tick labels
#ax2.yticks(fontsize=17)  # Y-axis tick labels
#plt.xscale("log")
ax21.set_xlabel(r'$g/\kappa_{L}$',fontsize = LABEL_FS)
#ax21.axvspan(0, 2.465, facecolor='b', alpha=0.5)
#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  

ax21.tick_params(labelsize=TICK_FS)  # font size of tick labels 
ax21.text(0.10, 0.89, '(b)', transform=ax21.transAxes, fontsize=PANEL_FS, fontweight='bold')

ax21.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left',
    bbox_to_anchor=(0.03, 0.23)
)

ax21.set_xscale("log")
#fig.supylabel("Cantidades termodinámicas", fontsize=22)
#plt.subplots_adjust(left=0.05) 
plt.tight_layout(pad=0.4)
plt.savefig("figthermog.pdf")
plt.show()
plt.close()

