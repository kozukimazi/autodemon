import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import cmath
import os


data = np.load("phononeVg=5_10^{-3}.npz")
eVs, entropf, Slr,cohes, concv, Nls,Nds,coheveig = data["eVs"], data["entropf"], data["Slr"], data["cohes"], data["concv"], data["Nls"], data["Nds"], data["coheveig"]
datasec = np.load("phononeVg=1.npz")
eVssec, entropfsec, Slrsec,cohessec, concvsec, Nlssec,Ndssec,coheveigsec,Qlrsec,Elrsec,Wtsec,Flrsec,Tislsec,Tidsec,Edssec,Fdsec,Wdfsec = datasec["eVs"], datasec["entropf"], datasec["Slr"], datasec["cohes"], datasec["concv"], datasec["Nls"], datasec["Nds"], datasec["coheveig"],datasec["Qlr"], datasec["Elr"],datasec["Wt"], datasec["Flr"], datasec["Tisl"], datasec["Tid"], datasec["Eds"], datasec["Fd"], datasec["Wdf"]

datasec2 = np.load("phononeVg=10.npz")
eVssec2, entropfsec2, Slrsec2,cohessec2, concvsec2, Nlssec2,Ndssec2,coheveigsec2,Qlrsec2,Elrsec2,Wtsec2,Flrsec2,Tislsec2,Tidsec2,Edssec2,Fdsec2,Wdfsec2 = datasec2["eVs"], datasec2["entropf"], datasec2["Slr"], datasec2["cohes"], datasec2["concv"], datasec2["Nls"], datasec2["Nds"], datasec2["coheveig"],datasec2["Qlr"], datasec2["Elr"],datasec2["Wt"], datasec2["Flr"], datasec2["Tisl"], datasec2["Tid"], datasec2["Eds"], datasec2["Fd"], datasec2["Wdf"]


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 



LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9


# Create subplots (1 row, 2 columns)
fig, (ax11, ax21) = plt.subplots(2, 1,sharex=True, figsize=(3.39, 4.5))  # 1 row, 2 columns

ax11.plot(eVs,Nls,color='black',lw=LINE_W, label = r'$g/\kappa_{L}= 5 \times 10^{-1}$')
ax11.plot(eVs,Nlssec,color='red',lw=LINE_W, label = r'$g/\kappa_{L}= 10^{2}$')
#ax11.plot(eVs,Nlssec2,color='red',lw=LINE_W, label = r'$g_{L}/\kappa_{L}= 10^{3}$')
ax11.plot(eVs,Nds, color='red',linestyle = '--',lw=LINE_W)
#(0.7,0.48)
#ax11.set_ylabel(r'$I/\kappa_{L}$',fontsize = 20)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
#ax10.xticks(fontsize=17)  
#ax10.yticks(fontsize=17)
#ax11.legend(bbox_to_anchor=(0.20, 0.98),fontsize=15,loc = "upper left", ncol = 2)
ax11.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax11.tick_params(labelbottom=False,direction='in', which='both',labelsize = TICK_FS)
ax11.text(0.07, 0.94, '(a)', transform=ax11.transAxes, fontsize=PANEL_FS, fontweight='bold')
ax11.set_ylabel(r'$\dot{N}_{L}/\kappa_{L}$',fontsize = LABEL_FS)

ax11.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left',
    bbox_to_anchor=(0.50, 0.30)
)

ax21.plot(eVs,coheveig,label = r'$\mathcal{C}_{l_{1}}(\mathrm{eig})(\mathrm{partial})$', color = 'black',linestyle = '--',lw = LINE_W)
ax21.plot(eVssec,coheveigsec,label = r'$\mathcal{C}_{l_{1}}(\mathrm{eig})(\mathrm{global})$', color = 'red',lw = LINE_W)    
#ax21.plot(eVssec2,coheveigsec2,label = r'$\mathcal{C}_{l_{1}}(\mathrm{eig})$', color = 'red',lw = LINE_W)
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
ax21.plot(eVs,concv, label = r'$\mathcal{C}_{on}(\mathrm{partial})$', color = 'black',lw=LINE_W) 
ax21.plot(eVssec,concvsec, label = r'$\mathcal{C}_{on}(\mathrm{global})$', color = 'green',linestyle = '--',lw=LINE_W)
#ax21.plot(eVssec2,concvsec2, label = r'$\mathcal{C}_{on}(\mathrm{sec})$', color = 'red',linestyle = '--',lw=LINE_W)
#plt.plot(eVs,Qdf,label = r'$J_{d}$',color = "gray",lw=2)
#ax2.xticks(fontsize=17)  # X-axis tick labels
#ax2.yticks(fontsize=17)  # Y-axis tick labels
#plt.xscale("log")
ax21.set_xlabel(r'$eV/T$',fontsize = LABEL_FS)
ax21.axvspan(0, 2.465, facecolor='b', alpha=0.5)

#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  

ax21.tick_params(labelsize=TICK_FS)  # font size of tick labels 
ax21.text(0.07, 0.94, '(b)', transform=ax21.transAxes, fontsize=PANEL_FS, fontweight='bold')

ax21.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left',
    bbox_to_anchor=(0.55, 0.95)
)

#fig.supylabel("Cantidades termodinámicas", fontsize=22)
#plt.subplots_adjust(left=0.05) 
plt.tight_layout(pad=0.4)
plt.savefig("figcurrent2.pdf")
plt.show()
plt.close()




#Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(3.39, 4.6))  # 1 row, 2 columns


LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

ax20.plot(eVs,Elrsec2,color='black',linestyle='--',lw=LINE_W, label = r'$\dot{E}_{2}$')
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
ax20.plot(eVs,Flrsec2,color='black',lw=LINE_W, label = r'$\dot{F}_{2}$')
ax20.plot(eVs,Tislsec2,label = r'$T\dot{I}_{2}$',linestyle='--', dashes=(5, 3), color = 'r',lw=LINE_W)
ax20.plot(eVs,Wtsec2,label = r'$\dot{W}_{2}$', color = 'm',lw=LINE_W)
ax20.plot(eVs,Qlrsec2, color='red',lw = LINE_W,label = r'$\dot{Q}_{2}$')
#ax10.xticks(fontsize=17)  
#ax10.yticks(fontsize=17)
ax20.set_xlabel(r'$eV/T$',fontsize = LABEL_FS)
ax20.tick_params(direction='in', which='both',labelsize=TICK_FS)

ax20.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax20.text(0.92, 0.1, '(b)', transform=ax20.transAxes, fontsize=PANEL_FS, fontweight='bold')

ax20.legend(
    fontsize=7,
    frameon=True,
    ncol=3,
    loc='upper left'
)



ax10.plot(eVs,Edssec, color='black',linestyle='--',lw=LINE_W, label = r'$\dot{E}_{1}$')
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='red',lw=2, label = r'$\dot{I}_{rl}$')
ax10.plot(eVs,Fdsec, color='black',lw=LINE_W, label = r'$\dot{F}_{1}$')
ax10.plot(eVs,Tidsec,label = r'$T_{B_D}\dot{I}_{1}$', linestyle='--', dashes=(5, 3),color = 'r',lw=LINE_W)
ax10.plot(eVs,Wdfsec,label = r'$\dot{W}_{1}$', color = 'm',lw=LINE_W)
#plt.plot(eVs,Qdf,label = r'$J_{d}$',color = "gray",lw=2)
#ax2.xticks(fontsize=17)  # X-axis tick labels
#ax2.yticks(fontsize=17)  # Y-axis tick labels
#plt.xscale("log")

#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  
ax10.tick_params(labelbottom=False,direction='in', which='both', labelsize=TICK_FS) # font size of tick labels 
ax10.text(0.92, 0.1, '(a)', transform=ax10.transAxes, fontsize=PANEL_FS, fontweight='bold')
ax10.axvspan(0, 2.465, facecolor='b', alpha=0.5)
#fig.supylabel("Cantidades termodinámicas", fontsize=22)
#plt.subplots_adjust(left=0.05) 

ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=2,
    loc='center left',
    bbox_to_anchor=(-0.005, 0.66)
)

ax10.text(
    0.50, 0.25,
    r'$g/\kappa_L=10^{2}$',
    transform=ax10.transAxes,
    fontsize=9.5
)



for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


plt.tight_layout(pad=0.4)
plt.savefig("figthermosec1.pdf")
plt.show()
plt.close()

