import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import cmath
import os
from scipy.sparse import kron, eye, csr_matrix


def anticonmutador(A,B):
    return np.matmul(A,B) + np.matmul(B,A)

def fermi(E,mu,beta):
    return 1/(np.exp((E-mu)*beta) + 1)


sigmax = np.array([[0,1],
                   [1,0]])

sigmay = np.array([[0,-1j],
                   [1j,0]])

iden = np.array([[1,0],
                 [0,1]])

sigmaz = np.array([[1,0],
                   [0,-1]])

sigmaup = (sigmax + 1j*sigmay)/2
sigmadown = (sigmax - 1j*sigmay)/2

auxdag = np.kron(sigmaz,sigmaup)
aux = np.kron(sigmaz,sigmadown)

auxd = np.kron(sigmaz,sigmaz)
#Jordan-Wigner
dldag = np.kron(sigmaup,np.eye(4))
dl = np.kron(sigmadown,np.eye(4))

drdag = np.kron(auxdag,np.eye(2))
dr = np.kron(aux,np.eye(2))

dddag = np.kron(auxd,sigmaup)
dd = np.kron(auxd,sigmadown)

tot = anticonmutador(dddag,dd)
totl = anticonmutador(dldag,dl)
totr = anticonmutador(drdag,dr)


nd = np.matmul(dddag,dd)
nl = np.matmul(dldag,dl)
nr = np.matmul(drdag,dr)

#print(nl)
#print(nd)
v0 = np.array([[0],
                 [1]])

v00 = np.array([[0],
                 [0],
                 [0],
                 [0],
                 [0],
                 [0],
                 [0],
                 [1]])

v100 = dldag @ v00
v010 = drdag @ v00
v001 = dddag @ v00
totrl = np.matmul(drdag,dldag)
totdl = np.matmul(dddag,dldag)
totdr = np.matmul(dddag,drdag)
totf = np.matmul(dddag,totrl )
v110 = totrl @ v00
v101 = totdl @ v00
v011 = totdr @ v00
v111 = totf @ v00

def partial_trace_d(A):
    """
    Partial trace over the third qubit (subsystem d).
    Input:  A — 8x8 matrix (operator on 3 qubits)
    Output: 4x4 matrix (operator on qubits l and r)
    """
    A = A.reshape(2, 2, 2, 2, 2, 2)  # reshape to (l, r, d, l', r', d')
    return np.trace(A, axis1=2, axis2=5).reshape(4,4)  # trace over d and d'

# d is MSB, then r, then l
def partial_trace_d_custom(A):
    A = A.reshape(2, 2, 2, 2, 2, 2)  # (d, r, l, d', r', l')
    # trace out d and d' → axes 0 and 3
    return np.trace(A, axis1=0, axis2=3).reshape(4,4)


alp,alp2,alp3,alp4,alp5 = 0.,0.,0.0,0.0,0.
a,b,c,d = 1j*alp,1j*alp2,1j*alp3,1j*alp4
#111,110,101,100,011,010,001,000
rho0 = np.array([[1/4,0,0,0,0,0,0,0],
                 [0,1/4,0,0,0,0,0,0],
                 [0,0,1/8,0,0.5,0,0,0],
                 [0,0,0,1/4,0,1,0,0],
                 [0,0,0.5,0,1,0,0,0],
                 [0,0,0,1,0,1/8,0,0],
                 [0,0,0,0,0,0,1,0],
                 [0,0,0,0,0,0,0,1/4]])


rhof = partial_trace_d(rho0)
print(rhof)

