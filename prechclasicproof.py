import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import linalg as la
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



#Jordan-Wigner
dsdag = np.kron(sigmaup,np.eye(2))
ds = np.kron(sigmadown,np.eye(2))

dddag = np.kron(sigmaz,sigmaup)
dd = np.kron(sigmaz,sigmadown)



##Fock Basis of three quantum dots
v00 = np.array([[0],
                [0],
                [0],
                [1]])

v10 = dsdag @ v00
v01 = dddag @ v00
totrl = np.matmul(dddag,dsdag)
v11 = totrl @ v00


basis = np.array([v00,v10,v01,v11])


nd = np.matmul(dddag,dd)
ns = np.matmul(dsdag,dd)
#perturbation superoperator

def pert(V,hbar=1):
    d = np.shape(V)[0]
    superH = (-1j/hbar) * (np.kron(np.eye(d), V ) - np.kron(V.T,  np.eye(d))  ) 
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
def Drazinalg(L0,tol):
    N = len(L0)
    #we need to order first such that the last eigenvalue is zero
    f1 = lambda x: abs(x) > tol 
    #The schur decomposition
    #A=QTQ^{-1}
    #T1 triangular matrix,Q1 the transformation
    #k1 specify the number of eigenvalues that f is true
    T1,Q1,k1 = la.schur(L0, sort=f1) 

    #Now we order the eigenvalues such that the first eigenvalue is zero
    f2 = lambda x: abs(x) < tol
    T2,Q2,k2 = la.schur(L0, sort=f1) 

    #The matrix transformation
    U = np.zeros((N,N),dtype=np.complex_)
    for i in range(k1):
        U[:,i] = Q1[:,i]

    for j in range(k1,N-k1):
        U[:,j] = Q2[:,j] 

    U1 = np.inv(U)
    V = np.natmul(U1,np.matmul(L0,U))
    #The descomposition here, the nonsigular matrix
    M = np.zeros((k1,k1),dtype=np.complex_)
    #The descomposition
    Z = np.zeros((N,N),dtype=np.complex_)
    #here we get M from the block of V
    if (k1!=0):
        for i in range(k1):
           for j in range(k1):               
                M[i,j] = V[i,j]  

    M1 = np.inv(M)
    #Here the Drazin
    for i in range(k1):
           for j in range(k1):               
                Z[i,j] = M1[i,j]                  

    Draz = np.matmul(U,np.matmul(Z,U1))
    return Draz

def Ds(E,U,mus,betas,gammas):
    Ns = np.matmul(dsdag,ds)
    Nd = np.matmul(dddag,dd)
    d = len(Ns)
    auxs1 = np.sqrt( fermi(E,mus,betas)*gammas )*np.matmul( (np.eye(d)-Nd),dsdag )
    auxs2 = np.sqrt( (1-fermi(E,mus,betas))*gammas )*np.matmul( (np.eye(d)-Nd),ds)
    auxs3 = np.sqrt( fermi(E+U,mus,betas)*gammas )*np.matmul( Nd,dsdag )
    auxs4 = np.sqrt( (1-fermi(E+U,mus,betas))*gammas )*np.matmul( Nd,ds)

    return [auxs1,auxs2,auxs3,auxs4]

def Dd(E,U,mud,betad,gammad):
    Ns = np.matmul(dsdag,ds)
    d = len(Ns)
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*np.matmul( (np.eye(d)-Ns),dddag )
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*np.matmul( (np.eye(d)-Ns),dd )
    auxd3 = np.sqrt( fermi(E+U,mud,betad)*gammad )*np.matmul( Ns,dddag )
    auxd4 = np.sqrt( (1-fermi(E+U,mud,betad))*gammad )*np.matmul( Ns,dd)   

    return [auxd1,auxd2,auxd3,auxd4]

def Ds(E,U,mus,betas,gammas):
    Ns = np.matmul(dsdag,ds)
    Nd = np.matmul(dddag,dd)
    d = len(Ns)
    auxs1 = np.sqrt( fermi(E,mus,betas)*gammas )*np.matmul( (np.eye(d)-Nd),dsdag )
    auxs2 = np.sqrt( (1-fermi(E,mus,betas))*gammas )*np.matmul( (np.eye(d)-Nd),ds)
    auxs3 = np.sqrt( fermi(E+U,mus,betas)*gammas )*np.matmul( Nd,dsdag )
    auxs4 = np.sqrt( (1-fermi(E+U,mus,betas))*gammas )*np.matmul( Nd,ds)

    return [auxs1,auxs2,auxs3,auxs4]

def Dd(E,U,mud,betad,gammad):
    Ns = np.matmul(dsdag,ds)
    d = len(Ns)
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*np.matmul( (np.eye(d)-Ns),dddag )
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*np.matmul( (np.eye(d)-Ns),dd )
    auxd3 = np.sqrt( fermi(E+U,mud,betad)*gammad )*np.matmul( Ns,dddag )
    auxd4 = np.sqrt( (1-fermi(E+U,mud,betad))*gammad )*np.matmul( Ns,dd)   

    return [auxd1,auxd2,auxd3,auxd4]


def Dissipator(E,U,mus,mud,betas,betad,gammas,gammad):
    DS = Ds(E,U,mus,betas,gammas)
    DD = Dd(E,U,mud,betad,gammad)

    tot = []
    for s in DS:
        tot.append(s)
    for d in DD:
        tot.append(d)    

    return tot

def Hamiltonian(E,U,g):
    a1 = E*ns + E*nd
    a2 = g*( np.matmul(dsdag,dd) + np.matmul(dddag,ds) )
    a3 = U*np.matmul(ns,nd) 
    return a1+a2+a3

def Inte(g):
    a2 = g*( np.matmul(dsdag,dd) + np.matmul(dddag,ds) )
    return a2

#revisar evolucion
def Evol(L0,L0draz,P,Q,V,rho0):
    #calcular PVQL^{-1}_{0}QVP
    paso1 = Prin(rho0)
    pasof0 = np.matmul(L0,paso1)
    paso2 = np.matmul(V,paso1)
    #revisar aqui
    paso3 = np.matmul(Q,paso2)
    paso4 = np.matmul(L0draz,paso3)
    paso5 = np.matmul(Q,paso4)
    paso6 = np.matmul(V,paso5)
    paso7 = Prin(paso6)
    Final = pasof0-paso7
    return 1

U0 = 0.
g0 = 0.005

eV = 6.5
mus1 = eV/2
mud1 = -eV/2
mul = eV/20

betas,betad,betal = 1,1,1
gs,gd,gl = 1/100,1/100,0
Ls = Dissipator(0,U0,mus1,mud1,betas,betad,gs,gd)
H = Hamiltonian(0,U0,g0)
V0 = Inte(g0)
Vnu = pert(V0)


rho0 = np.array([[1/4,0,0,0],
                 [0,1/4,0.1j,0],
                 [0,-0.1j,1/4,0],
                 [0,0,0,1/4]])



gss = np.linspace(0.,1,500)
gaux = []
Il = []
I2l = []
Nss = 3
for g in gss:
    t=4000
    gau = g/gs
    H0 = Hamiltonian(0,U0,g)
    gaux.append(gau)
    
    #print(Il0/gs)  
    #print(g)
    

plt.plot( gaux,Il)
plt.ylabel(r'$I_{L}/\gamma$',fontsize = 20)     
plt.xlabel(r'$g/\gamma$',fontsize = 20)
plt.xscale("log")
plt.show()
plt.plot( gaux,I2l)
plt.ylabel(r'$\langle \langle I^{2}_{L} \rangle \rangle/\gamma$',fontsize = 20)     
plt.xlabel(r'$g/\gamma$',fontsize = 20)
plt.xscale("log")
plt.show()

