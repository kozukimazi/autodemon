import numpy as np
import matplotlib.pyplot as plt

# Shared bath properties:
ep0 = 1/800
alp0 = 5/600
gamma = ep0 + alp0   # coupling strength
Uf = 40
W = np.sqrt((ep0/alp0)*(Uf**2))  # cut-off
T = 100  # temperature
beta = 1. / T

# Chemical potentials for the two baths:
mu_L = 10
mu_R = -10

E = 0
U = 40
omega0 = E+U

def Lorent(gamf,W,e0,w):

    tot = gamf*W**2
    totf = (w-e0)**2 + W**2
    return tot/totf

lore = []
omegas = np.linspace(-400,400,500)

for ws in omegas:

    lor = Lorent(gamma,W,E,ws)
    lore.append(lor)

plt.plot(omegas,lore)
plt.show()