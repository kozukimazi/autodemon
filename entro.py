import numpy as np
import matplotlib.pyplot as plt

data = np.load("entropy.npz")
eVs, entropf,Slr,Nds = data["eVs"], data["entropf"], data["Slr"], data["Nds"]

dataph = np.load("entropyph.npz")
eVsph, entropfph,Slrt,Ndsph = dataph["eVs"], dataph["entropf"], dataph["Slrt"], dataph["Nds"]

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 

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
ax10.plot(eVs, entropf, color='black',   lw=LINE_W, label=r'$\dot{\sigma}^{O}_{LR}$')
#ax10.plot(Jof1, Imalphg1, color='orange', lw=LINE_W, label=r'$10^{-2}$')
ax10.plot(eVs, Slr, color='red',  lw=LINE_W, label=r'$\dot{\sigma}_{LR}$')
ax10.plot(eVs, Nds, color='black',lw=LINE_W,linestyle = '--')
#ax10.set_ylabel(r'$2g\mathrm{Im}(\beta)$', fontsize=LABEL_FS)
#ax10.set_xlabel(r'$eV/T$', fontsize=LABEL_FS)
#ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax10.text(0.9, 0.05, '(a)', transform=ax10.transAxes,
          fontsize=PANEL_FS, fontweight='bold')
ax10.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax10.legend(
    fontsize=7,
    frameon=True,
    ncol=1,
    loc='upper left'
)

# ---------- Panel (b)
ax20.plot(eVsph, entropfph, color='black',   lw=LINE_W)
#ax20.plot(Jof1, Imbetg1, color='orange', lw=LINE_W)
ax20.plot(eVsph, Slrt, color='red',  lw=LINE_W)
ax20.plot(eVsph, Ndsph, color='black',lw=LINE_W, linestyle = '--')
ax20.set_xlabel(r'$eV/T$', fontsize=LABEL_FS)
#ax20.set_ylabel(r'$2g\,\mathrm{Im}(\alpha)$', fontsize=LABEL_FS)
#ax20.set_xscale('log')
ax20.tick_params(direction='in', which='both', labelsize=TICK_FS)
ax20.text(0.9, 0.05, '(b)', transform=ax20.transAxes,
          fontsize=PANEL_FS, fontweight='bold')
ax20.axvspan(0, 2.465, facecolor='b', alpha=0.5)
ax20.text(
    0.05, 0.30,
    r'$J_0/(\beta_{\mathrm{ph}}\kappa_L)=10^{-1}$',
    transform=ax20.transAxes,
    fontsize=8
)

# ---------- Spines
for ax in (ax10, ax20):
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)

plt.savefig("entropyphtot.pdf")
plt.close()