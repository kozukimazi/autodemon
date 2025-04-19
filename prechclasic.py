import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import cmath
import os


def anticonmutador(A,B):
    return np.matmul(A,B) + np.matmul(B,A)

def fermi(E,mu,beta):
    return 1/(np.exp((E-mu)*beta) + 1)

def derivada(t,x):
    N = np.shape(t)[0]
    der = []
    ts = []
    for i in range(N-1):
        derivada = (x[i+1] - x[i] )/(t[i+1] - t[i])
        der.append(derivada)
        ts.append(t[i])

    return ts,der  

def quadrature(x1,y1):
    n = len(x1)-1
    total = 0
    for ns in range(n):
        total += (x1[ns+1] - x1[ns])*(y1[ns+1] + y1[ns])*(1/2)
    return total  

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

##Fock Basis of three quantum dots
v00 = np.array([[0],
                [0],
                [0],
                [0],
                [0],
                [0],
                [0],
                [1]])

v100 = dldag @ v00
v010 = -drdag @ v00
v001 = dddag @ v00
totrl = np.matmul(drdag,dldag)
totdl = np.matmul(dddag,dldag)
totdr = np.matmul(dddag,drdag)
totf = np.matmul(dddag,totrl )
v110 = totrl @ v00
v101 = -totdl @ v00
v011 = totdr @ v00
v111 = totf @ v00

basis = np.array([v00,v100,v010,v001,v110,v101,v011,v111])

#perturbation superoperator

def pert(g,V,hbar=1):
    d = np.shape(V)[0]
    superH = g*(-1j/hbar) * (np.kron(np.eye(d), V ) - np.kron(V.T,  np.eye(d))  ) 
    return superH


#we want to vectorize the P operator so we define
#We first vectorize rho0
def Vectorization(rho0):
    d = len(rho0)
    vec_rho =  np.reshape(rho0,(d**2,1))
    return vec_rho

#Here we vectorize the principal Part, rho0 is the vectorization
def Prin(rho0,basis):
    d=len(rho0)
    P = np.zeros(d**2,1)
    n = len(basis)
    for i in range(n):
        coef = (np.kron(basis[i].conj(),basis[i]).T)@rho0
        P += np.matmul(coef,np.kron(basis[i].conj(),basis[i]))
    return P    

##The coherence part

def Qpart(rho0,Princ):
    return rho0-Princ

#The liouvillian part

def Liouvillian( H,Ls, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    superL = sum( [np.kron(L.conjugate(),L) - 1/2 * (np.kron( np.eye(d), L.conjugate().T.dot(L)) +
                                                     np.kron( L.T.dot(L.conjugate()),np.eye(d) ))
                                                      for L in Ls ] )    
    return superH + superL

#Here we need to calculate the Drazin part

