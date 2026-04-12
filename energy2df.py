import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import matplotlib.colors as colors

data = np.load("energy.npz")
x = data["x"]
y = data["y"]
Z = data["Z"]

data2 = np.load("energyfs.npz")
x2 = data2["x"]
y2 = data2["y"]
Z2 = data2["Z"]



plt.imshow(Z, cmap='viridis')
plt.colorbar()
plt.show()

plt.imshow(Z2, cmap='viridis')
plt.colorbar()
plt.show()



Z_total = np.vstack((Z, Z2))
y_total = np.concatenate((y, y2))

plt.imshow(Z_total, cmap='viridis')
plt.colorbar()
plt.show()


plt.imshow(
    Z_total,
    extent=[min(x), max(x), min(y_total), max(y_total)],
    origin='lower',
    aspect='auto'
)
plt.colorbar()
plt.xlabel(r'$J_0/(\beta_{\mathrm{ph}} \kappa_L)$')
plt.ylabel(r'$g/ \kappa_L$')
plt.show()


###mejorar formato############
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif" 

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300
})


fig, ax = plt.subplots(figsize=(3.4, 3))

pcm = ax.pcolormesh(x, y_total, Z_total, shading='auto', cmap='viridis')

cbar = plt.colorbar(pcm, ax=ax)
cbar.set_label(r'$\dot{E}_1$')

ax.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}} \kappa_L)$')
ax.set_ylabel(r'$g/ \kappa_L$')

# 🔥 aquí está la clave
ax.set_xscale("log")
ax.set_yscale("log")

plt.tight_layout()
plt.show()


import matplotlib.ticker as ticker
import matplotlib.colors as colors



plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300
})

fig, ax = plt.subplots(figsize=(3.4, 3))

pcm = ax.pcolormesh(
    x, y_total, Z_total,
    shading='auto',
    cmap='viridis'
)

cbar = plt.colorbar(pcm, ax=ax)
cbar.set_label(r'$\dot{E}_1$')
cbar.ax.tick_params(labelsize=8, direction='in', width=0.8)
cbar.outline.set_linewidth(0.8)

ax.set_xlabel(r'$J_0/(\beta_{\mathrm{ph}} \kappa_L)$')
ax.set_ylabel(r'$g/ \kappa_L$')

ax.set_xscale("log")
ax.set_yscale("log")

ax.xaxis.set_major_locator(ticker.LogLocator(base=10))
ax.yaxis.set_major_locator(ticker.LogLocator(base=10))

ax.tick_params(which='both', direction='in', width=0.8)

for spine in ax.spines.values():
    spine.set_linewidth(0.8)

plt.tight_layout(pad=0.3)
plt.savefig("figure3D.pdf", bbox_inches='tight')
plt.show()