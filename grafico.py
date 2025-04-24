import matplotlib.pyplot as plt
import numpy as np

def lorentzian(e,w,gamma):

    val = (e-w)**2 + gamma**2
    return 1/val

num = 5000
eps = np.linspace(0,100,num)

g1 = []
g2 = []
gaux = []
gaux2 = []
U = 50
epsilon = 25

for es in eps:
    v1 = lorentzian(es,epsilon + U,5)
    v2 =lorentzian(es,epsilon,5)
    c = 0.02
    g1.append(v1+c)
    g2.append(v2+c)
    gaux.append(0)
    gaux2.append(c)


plt.plot(eps,g1,color='red',lw = 4, label = r'$\gamma_{R}(\epsilon)$')
plt.plot(eps,g2,color='blue',lw = 4, label = r'$\gamma_{L}(\epsilon)$')
plt.plot(eps,gaux,color='black',lw = 4)
plt.plot(eps,gaux2,linestyle='--', dashes=(5, 9),color='black',lw = 4)
plt.xticks([])  # X-axis tick labels
plt.yticks([])
plt.legend(fontsize = 15)
plt.show()    

