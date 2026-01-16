import os
import matplotlib.pyplot as plt
import numpy as np

g_data = []
coheparcial = []
concuparcial = []
Qlrparcial = []
Wlrparcial = []
Ilrparcial = []

cohered = []
concured = []
Qlrred = []
Wlrred = []

coheglobal = []
concuglobal = []
Qlrglobal = []
Wlrglobal = []
Ilrglobal = []

#######################################
########compararcodigos################
#datos = 'parcialweakg'
datos = 'semicouplarge'
#datos = 'parcialweakg'
fichero = open(datos)
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    g_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    concuparcial.append(float(item)) 
    #print(float(item))  
nuevo.close()

nuevo = open(datos)
for item in [data.split()[2] for data in nuevo]:  
    coheparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()

nuevo = open(datos)
for item in [data.split()[3] for data in nuevo]:  
    Qlrparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()

nuevo = open(datos)
for item in [data.split()[4] for data in nuevo]:  
    Wlrparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()

nuevo = open(datos)
for item in [data.split()[5] for data in nuevo]:  
    Ilrparcial.append(float(item))  
    #print(float(item)) 
nuevo.close()

#######################################
########compararcodigos################
#datos = 'parcialweakg'
datos0 = 'redcouplarge'
#datos = 'parcialweakg'


nuevo0 = open(datos0)
for item in [data.split()[1] for data in nuevo0]:  
    concured.append(float(item)) 
    #print(float(item))  
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[2] for data in nuevo0]:  
    cohered.append(float(item))  
    #print(float(item)) 
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[3] for data in nuevo0]:  
    Qlrred.append(float(item))  
    #print(float(item)) 
nuevo0.close()

nuevo0 = open(datos0)
for item in [data.split()[4] for data in nuevo0]:  
    Wlrred.append(float(item))  
    #print(float(item)) 
nuevo0.close()


#datos = 'parcialweakg'
datos1 = 'globalcouplarge'
#datos = 'parcialweakg'


nuevo1 = open(datos1)
for item in [data.split()[1] for data in nuevo1]:  
    concuglobal.append(float(item)) 
    #print(float(item))  
nuevo1.close()

nuevo1 = open(datos1)
for item in [data.split()[2] for data in nuevo1]:  
    coheglobal.append(float(item))  
    #print(float(item)) 
nuevo1.close()

nuevo1 = open(datos1)
for item in [data.split()[3] for data in nuevo1]:  
    Qlrglobal.append(float(item)) 
    #print(float(item))  
nuevo1.close()

nuevo1 = open(datos1)
for item in [data.split()[4] for data in nuevo1]:  
    Wlrglobal.append(float(item))  
    #print(float(item)) 
nuevo1.close()

nuevo1 = open(datos1)
for item in [data.split()[4] for data in nuevo1]:  
    Ilrglobal.append(float(item))  
    #print(float(item)) 
nuevo1.close()

plt.plot(g_data,coheparcial,label = "parcial", color = 'b',lw = 2)
plt.plot(g_data,cohered, linestyle='--', label = "redfield", color = 'r',lw=2)  
plt.plot(g_data,coheglobal, linestyle='--', label = "global", color = 'black',lw=2)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel("Coherencia",fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()


plt.plot(g_data,concuparcial,label = "parcial", color = 'b',lw = 2)
plt.plot(g_data,concured, linestyle='--', label = "redfield", color = 'r',lw=2)  
plt.plot(g_data,concuglobal, linestyle='--', label = "global", color = 'black',lw=2)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel("Concurrencia",fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()

plt.plot(g_data,Qlrparcial,label = "parcial", color = 'b',lw = 3)
plt.plot(g_data,Qlrglobal, linestyle='--', label = "global", color = 'black',lw=3)  
plt.plot(g_data,Qlrred,linestyle='--', label = "redfield", color = 'red',lw=3)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel(r'$J_{LR}$',fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()

plt.plot(g_data,Wlrparcial,label = "parcial", color = 'b',lw = 3)
plt.plot(g_data,Wlrglobal, linestyle='--', label = "global", color = 'black',lw=3)  
plt.plot(g_data,Wlrred,linestyle='--', label = "redfield", color = 'red',lw=3)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel(r'$W_{LR}$',fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()

plt.plot(g_data,Ilrparcial,label = "parcial", color = 'b',lw = 3)
plt.plot(g_data,Ilrglobal, linestyle='--', label = "global", color = 'black',lw=3)  
plt.xlabel(r'$g/\kappa_{L}$',fontsize = 22)   
plt.ylabel(r'$\dot{I}_{LR}$',fontsize = 22) 
plt.xscale("log")
plt.xticks(fontsize=22)  
plt.yticks(fontsize=22)
plt.legend(fontsize=19,loc = "upper left")
plt.show()


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 


data0 = np.load("phononJ0=5_10^{-1}comp.npz")
gof10,ql0,qph0,Id0,Ile0,Ire0,Iphs0,cohes0,concv0,Nls0 = data0["gof1"], data0["Qlr"],data0["Qphs"], data0["Id"], data0["Ile"], data0["Ire"], data0["Iphs"], data0["cohes"], data0["concv"], data0["Nls"]

data1 = np.load("phononJ0=5_10^{-1}redcomp.npz")
gof11,ql1,qph1,cohes1,concv1,Nls1 = data1["gof1"], data1["Qlrs"], data1["Qphs"], data1["cohev"], data1["concuv"], data1["Nls"]

data2 = np.load("phononJ0=5_10^{-1}seccomp.npz")
gof12,ql2,qph2,cohes2,concv2,Nls2 = data2["gof1"], data2["Qlrs"], data2["Qphs"], data2["cohev"], data2["concuv"], data2["Nls"]


data0x = np.load("phononJ0=10comp.npz")
gof10x,ql0x,qph0x,Id0x,Ile0x,Ire0x,Iphs0x,cohes0x,concv0x,Nls0x = data0x["gof1"], data0x["Qlr"],data0x["Qphs"], data0x["Id"], data0x["Ile"], data0x["Ire"], data0x["Iphs"], data0x["cohev"], data0x["concv"], data0x["Nls"]

data1x = np.load("phononJ0=10redcomp.npz")
gof11x,ql1x,qph1x,cohes1x,concv1x,Nls1x = data1x["gof1"], data1x["Qlrs"], data1x["Qphs"], data1x["cohev"], data1x["concuv"], data1x["Nls"]

data2x = np.load("phononJ0=10seccomp.npz")
gof12x,ql2x,qph2x,cohes2x,concv2x,Nls2x = data2x["gof1"], data2x["Qlrs"], data2x["Qphs"], data2x["cohev"], data2x["concuv"], data2x["Nls"]


data0y = np.load("phononJ0=10^{-3}comp.npz")
gof10y,ql0y,qph0y,Id0y,Ile0y,Ire0y,Iphs0y,cohes0y,concv0y,Nls0y = data0y["gof1"], data0y["Qlr"],data0y["Qphs"], data0y["Id"], data0y["Ile"], data0y["Ire"], data0y["Iphs"], data0y["cohev"], data0y["concv"], data0y["Nls"]

data1y = np.load("phononJ0=10^{-3}redcomp.npz")
gof11y,ql1y,qph1y,cohes1y,concv1y,Nls1y = data1y["gof1"], data1y["Qlrs"], data1y["Qphs"], data1y["cohev"], data1y["concuv"], data1y["Nls"]

data2y = np.load("phononJ0=10^{-3}seccomp.npz")
gof12y,ql2y,qph2y,cohes2y,concv2y,Nls2y = data2y["gof1"], data2y["Qlrs"], data2y["Qphs"], data2y["cohev"], data2y["concuv"], data2y["Nls"]

n0f = len(Nls2y) 
gl = 1/100
Nls0f = []
Nls1f = []
Nls2f = []
Nls0fx = []
Nls1fx = []
Nls2fx = []
Nls0fy = []
Nls1fy = []
Nls2fy = []
for i in range(n0f):
    Nls0f.append(Nls0[i]/gl)
    Nls1f.append(Nls1[i]/gl)
    Nls2f.append(Nls2[i]/gl)
    Nls0fx.append(Nls0x[i]/gl)
    Nls1fx.append(Nls1x[i]/gl)
    Nls2fx.append(Nls2x[i]/gl)
    Nls0fy.append(Nls0y[i]/gl)
    Nls1fy.append(Nls1y[i]/gl)
    Nls2fy.append(Nls2y[i]/gl)


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
ax10.plot(gof10, cohes0, color='blue',   lw=LINE_W, label="semilocal")
#ax10.plot(Jof1, Imalphg1, color='orange', lw=LINE_W, label=r'$10^{-2}$')
ax10.plot(gof12, cohes2, color='red',    lw=LINE_W, label="global")
ax10.plot(gof11, cohes1, color='green',  lw=LINE_W, label="redfield", ls = '--')
ax10.plot(gof10x, cohes0x, color='blue',   lw=LINE_W, label="semilocal", ls = '--')
ax10.plot(gof12x, cohes2x, color='red',    lw=LINE_W, label="global",ls = '--')
ax10.plot(gof11x, cohes1x, color='green',  lw=LINE_W, label="redfield", ls = ':')

ax10.plot(gof10y, cohes0y, color='blue',   lw=LINE_W, label="semilocal", ls = '-.')
ax10.plot(gof12y, cohes2y, color='red',    lw=LINE_W, label="global",ls = '-.')
ax10.plot(gof11y, cohes1y, color='green',  lw=LINE_W, label="redfield", ls = '-.')


ax10.set_ylabel(r'$\mathcal{C}_{l_1}$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax10.text(0.05, 0.97, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='center left'
)

# ---------- Panel (b)
ax20.plot(gof10, Nls0, color='blue',   lw=LINE_W)
#ax20.plot(Jof1, Imbetg1, color='orange', lw=LINE_W)
ax20.plot(gof12, Nls2, color='red',    lw=LINE_W)
ax20.plot(gof11, Nls1, color='green',  lw=LINE_W,ls = '--')

ax20.plot(gof10x, Nls0x, color='blue',   lw=LINE_W, ls = '--')
#ax20.plot(Jof1, Imbetg1, color='orange', lw=LINE_W)
ax20.plot(gof12x, Nls2x, color='red',    lw=LINE_W, ls = '--')
ax20.plot(gof11x, Nls1x, color='green',  lw=LINE_W,ls = ':')


ax20.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$\dot{N}_L$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.05, 0.97, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("figcomparacion.pdf")
plt.close()


#############################
##########J0/bpkl=10##############


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
ax10.plot(gof10x, cohes0x, color='blue',   lw=LINE_W, label="semilocal")
#ax10.plot(Jof1, Imalphg1, color='orange', lw=LINE_W, label=r'$10^{-2}$')
ax10.plot(gof12x, cohes2x, color='red',    lw=LINE_W, label="global")
ax10.plot(gof11x, cohes1x, color='green',  lw=LINE_W, label="redfield", ls = '--')

ax10.set_ylabel(r'$\mathcal{C}_{l_1}$', fontsize=LABEL_FS)
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
ax20.plot(gof10x, Nls0x, color='blue',   lw=LINE_W)
#ax20.plot(Jof1, Imbetg1, color='orange', lw=LINE_W)
ax20.plot(gof12x, Nls2x, color='red',    lw=LINE_W)
ax20.plot(gof11x, Nls1x, color='green',  lw=LINE_W,ls = '--')


ax20.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$\dot{N}_L$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("figcomparacion1.pdf")
plt.close()


fig, (ax10, ax20) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(3.39, 4.0)
)

# ---------- Panel (a)
ax10.plot(gof10y, cohes0y, color='blue',   lw=LINE_W, label="semilocal")
#ax10.plot(Jof1, Imalphg1, color='orange', lw=LINE_W, label=r'$10^{-2}$')
ax10.plot(gof12y, cohes2y, color='red',    lw=LINE_W, label="global")
ax10.plot(gof11y, cohes1y, color='green',  lw=LINE_W, label="redfield", ls = '--')

ax10.set_ylabel(r'$\mathcal{C}_{l_1}$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax10.text(0.91, 0.88, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

ax10.text(
    0.05, 0.95,
    r'$J_0/(\beta_{Ph} \kappa_L)= 10^{-3}$',
    transform=ax.transAxes,
    fontsize=9
)

ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left'
)

# ---------- Panel (b)
ax20.plot(gof10y, Nls0fy, color='blue',   lw=LINE_W)
#ax20.plot(Jof1, Imbet1, color='orange', lw=LINE_W)
ax20.plot(gof12y, Nls2fy, color='red',    lw=LINE_W)
ax20.plot(gof11y, Nls1fy, color='green',  lw=LINE_W,ls = '--')


ax20.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$\dot{N}_L/\kappa_L$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.91, 0.88, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.savefig("figcomparacion2.pdf")
plt.close()





# ---------------- Global style ----------------
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# ---------------- Figure geometry ----------------
fig, axs = plt.subplots(
    2, 2,
    figsize=(7.0, 4.4),   # APS 2-column width
    sharex='col'
)

(ax10, ax11), (ax20, ax21) = axs

LINE_W   = 2.2
LABEL_FS = 9
TICK_FS  = 8
PANEL_FS = 9
LEG_FS   = 7

# ===================== Panel (a) =====================
ax10.plot(gof10y, cohes0y, color='blue',  lw=LINE_W, label="semilocal")
ax10.plot(gof12y, cohes2y, color='red',   lw=LINE_W, label="global")
ax10.plot(gof11y, cohes1y, color='green', lw=LINE_W, ls='--', label="redfield")

ax10.set_ylabel(r'$\mathcal{C}_{l_1}$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)

ax10.text(0.92, 0.88, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

ax10.text(
    0.05, 0.65,
    r'$J_0/(\beta_{Ph}\kappa_L)=10^{-3}$',
    transform=ax10.transAxes,
    fontsize=8
)

ax10.legend(
    fontsize=LEG_FS,
    frameon=True,
    loc='upper left'
)

# ===================== Panel (b) =====================
ax20.plot(gof10y, Nls0fy, color='blue',  lw=LINE_W)
ax20.plot(gof12y, Nls2fy, color='red',   lw=LINE_W)
ax20.plot(gof11y, Nls1fy, color='green', lw=LINE_W, ls='--')

ax20.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax20.set_ylabel(r'$\dot{N}_L/\kappa_L$', fontsize=LABEL_FS)
ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)

ax20.text(0.92, 0.88, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')

# ===================== Panel (c) — placeholder =====================
# Add your third observable here
ax11.plot(gof10, cohes0, color='blue',  lw=LINE_W)
ax11.plot(gof12, cohes2, color='red',   lw=LINE_W)
ax11.plot(gof11, cohes1, color='green', lw=LINE_W, ls='--')
ax11.set_xscale('log')
ax11.text(0.92, 0.88, '(c)', transform=ax11.transAxes,
          fontsize=PANEL_FS, fontweight='bold')
ax11.tick_params(direction='in', which='both', labelsize=TICK_FS)

ax11.text(
    0.05, 0.65,
    r'$J_0/(\beta_{Ph}\kappa_L)=5 \times 10^{-1}$',
    transform=ax11.transAxes,
    fontsize=8
)


# ===================== Panel (d) — placeholder =====================
# Add your fourth observable here
ax21.plot(gof10, Nls0f, color='blue',  lw=LINE_W)
ax21.plot(gof12, Nls2f, color='red',   lw=LINE_W)
ax21.plot(gof11, Nls1f, color='green', lw=LINE_W, ls='--')
ax21.text(0.92, 0.88, '(d)', transform=ax21.transAxes,
          fontsize=PANEL_FS, fontweight='bold')
ax21.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax21.tick_params(direction='in', which='both', labelsize=TICK_FS)

# ---------------- Shared cosmetics ----------------
for ax in axs.flat:
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

# ---------------- Layout & save ----------------
plt.tight_layout(pad=0.6)
plt.savefig("fig_2x2_PR.pdf")
plt.close()
