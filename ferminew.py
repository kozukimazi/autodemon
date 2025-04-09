# Imports
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigvalsh



# Convenience functions and parameters:
def deltafun(j, k):
    # Kronecker delta function. 
    return 1.0 if j == k else 0.

def f_approx(x, Nk):
    # Padé approximation to Fermi distribution. #
    f = 0.5
    for ll in range(1, Nk + 1):
        # kappa and epsilon are calculated further down
        f = f - 2 * kappa[ll] * x / (x**2 + epsilon[ll]**2)
    return f

def kappa_epsilon(Nk):
    # Calculate kappa and epsilon coefficients. 

    alpha = np.zeros((2 * Nk, 2 * Nk))
    for j in range(2 * Nk):
        for k in range(2 * Nk):
            alpha[j][k] = (
                (deltafun(j, k + 1) + deltafun(j, k - 1))
                / np.sqrt((2 * (j + 1) - 1) * (2 * (k + 1) - 1))
            )

    eps = [-2. / val for val in eigvalsh(alpha)[:Nk]]

    alpha_p = np.zeros((2 * Nk - 1, 2 * Nk - 1))
    for j in range(2 * Nk - 1):
        for k in range(2 * Nk - 1):
            alpha_p[j][k] = (
                (deltafun(j, k + 1) + deltafun(j, k - 1))
                / np.sqrt((2 * (j + 1) + 1) * (2 * (k + 1) + 1))
            )

    chi = [-2. / val for val in eigvalsh(alpha_p)[:Nk - 1]]

    eta_list = [
        0.5 * Nk * (2 * (Nk + 1) - 1) * (
            np.prod([chi[k]**2 - eps[j]**2 for k in range(Nk - 1)]) /
            np.prod([
                eps[k]**2 - eps[j]**2 + deltafun(j, k) for k in range(Nk)
            ])
        )
        for j in range(Nk)
    ]

    kappa = [0] + eta_list
    epsilon = [0] + eps

    return kappa, epsilon

# Number of expansion terms to retain:
Nk = 12
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

kappa, epsilon = kappa_epsilon(Nk)

# Phew, we made it to function that calculates the coefficients for the
# correlation function expansions:

def C(sigma, mu,omega0, Nk):
# Calculate the expansion coefficients for C_\sigma.
    beta = 1. / T
    ck = [0.5 * gamma * W * f_approx(1.0j * beta * W + sigma * beta*(omega0-mu), Nk)]
    vk = [W - sigma * 1.0j * omega0]
    for ll in range(1, Nk + 1):
        ck.append(
            -1.0j * (kappa[ll] / beta) * gamma * W**2
            / (((1j*sigma*epsilon[ll]) / beta + (mu-omega0) )**2 + W**2)
        )
        vk.append(epsilon[ll] / beta - sigma * 1.0j * mu)
    return ck, vk

def Ct(Nk,ck,vk,t):
    tot = 0
    for i in range(Nk+1):
        tot += ck[i]*np.exp(-vk[i]*t)
    return tot.real, tot.imag    

E = 4
U = 40
omega0 = E+U

ck_plus_L, vk_plus_L = C(1.0, mu_L,E, Nk)  # C_+, left bath
ck_minus_L, vk_minus_L = C(-1.0, mu_L,E, Nk)  # C_-, left bath

ck_plus_R, vk_plus_R = C(1.0, mu_R,E+U, Nk)  # C_+, right bath
ck_minus_R, vk_minus_R = C(-1.0, mu_R,E+U, Nk)  # C_-, right bath

#Nk0 = 8
#ck_plus_L0, vk_plus_L0 = C(1.0, mu_L, Nk0)  # C_+, left bath
#ck_minus_L0, vk_minus_L0 = C(-1.0, mu_L, Nk0)  # C_-, left bath

#ck_plus_R0, vk_plus_R0 = C(1.0, mu_R, Nk0)  # C_+, right bath
#ck_minus_R0, vk_minus_R0 = C(-1.0, mu_R, Nk0)  # C_-, right bath

ts = np.linspace(0,0.7,1000)
csre = []
csimag = []
csreR = []
csimagR = []
for t in ts:
    #cre,cim = Ct(Nk,ck_plus_L,vk_plus_L,t)
    #creR,cimR = Ct(Nk,ck_plus_R,vk_plus_R,t)
    cre,cim = Ct(Nk,ck_minus_L,vk_minus_L,t)
    creR,cimR = Ct(Nk,ck_minus_R,vk_minus_R,t)
    #cre0,cim0 = Ct(Nk0,ck_plus_L0,vk_plus_L0,t)
    csre.append(cre)
    csimag.append(cim)
    csreR.append(creR)
    csimagR.append(cimR)
    #csre0.append(cre0)
    #csimag0.append(cim0)

plt.plot(ts,csre, color='blue',lw=3,label = r'$\text{Re}[C^{-}_{fL}(t)]$')
plt.plot(ts,csreR, color='red',lw=3,label = r'$\text{Re}[C^{-}_{fR}(t)]$')
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.xlabel(r'$t$',fontsize=25)
plt.legend(fontsize=15,loc = "upper left")
#plt.xscale("log")
plt.show()   

plt.plot(ts,csimag, color='blue',lw=3,label = r'$\text{Im}[C^{-}_{fL}(t)]$')
plt.plot(ts,csimagR, color='red',lw=3,label = r'$\text{Im}[C^{-}_{fR}(t)]$')
#plt.plot(ts,csimag0)
plt.xlabel(r'$t$',fontsize=25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper left")
#plt.xscale("log")
plt.show()   
