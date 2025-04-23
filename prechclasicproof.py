import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy.linalg import inv
from scipy.linalg import eig
from scipy import linalg as la
from scipy import integrate
import cmath
import os


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
    d = len(V)
    superH = (-1j/hbar) * (np.kron(np.eye(d), V ) - np.kron(V.T,  np.eye(d))  ) 
    return superH


#we want to vectorize the P operator so we define
#We first vectorize rho0
def Vectorization(rho0):
    d = len(rho0)
    vec_rho =  np.reshape(rho0,(d**2,1))
    return vec_rho

#vectorization back to matrix
def vecback(rhovec0,d):
    
    rhoff = np.reshape(rhovec0,(d,d))
    return rhoff

#Here we vectorize the principal Part, rho0 is the vectorization
def Prin(rho0,basis):
    d=len(rho0)
    n = len(basis)
    P = np.zeros((d**2,d**2),dtype=np.complex_)
    for i in range(n):
        #print((basis[i].T))
        coef = (np.kron(basis[i].conj(),basis[i]))
        coef2 = (np.kron(basis[i].T,basis[i].conj().T))
        P += np.matmul(coef,coef2)
    return P    

##The coherence part

def Qpart(rho0,P):
    d = len(rho0)
    Iden = np.eye(d**2)
    return Iden-P


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
    T2,Q2,k2 = la.schur(L0, sort=f2) 

    #The matrix transformation
    U = np.zeros((N,N),dtype=np.complex_)
    for i in range(k1):
        U[:,i] = Q1[:,i]

    for j in range(k1,N):
        U[:,j] = Q2[:,j-k1] 

    #print(U)
    U1 = inv(U)
    V = np.matmul(U1,np.matmul(L0,U))
    #The descomposition here, the nonsigular matrix
    M = np.zeros((k1,k1),dtype=np.complex_)
    #The descomposition
    Z = np.zeros((N,N),dtype=np.complex_)
    #here we get M from the block of V
    if (k1!=0):
        for i in range(k1):
           for j in range(k1):               
                M[i,j] = V[i,j]  

    M1 = inv(M)
    #Here the Drazin
    for i in range(k1):
           for j in range(k1):               
                Z[i,j] = M1[i,j]                  

    Draz = np.matmul(U,np.matmul(Z,U1))
    return Draz

###############################################
#########another form to do de Drazin##########
###############################################
def Drazinspectral(L0,tol):
    N = len(L0)
    w, eigenl, eigenr = eig(L0, left=True)
    #print(np.shape(eigenl))
    #print(np.shape(eigenr))
    
    L0_D = np.zeros((N,N),dtype = np.complex_) 
    #print(np.shape(L0_D))
    for i in range(len(w)):
        #print(w[i])
        if (abs(w[i])>tol):
            #print(type(w[i]))            
            L0_D += np.outer(eigenr[:, i], eigenl[:, i])/w[i]

    return L0_D


##########################################
##############parteespecificadelsistema###
##########################################
##########################################

def Ds(E,mus,betas,gammas):
    Ns = np.matmul(dsdag,ds)
    auxs1 = np.sqrt( fermi(E,mus,betas)*gammas )*dsdag 
    auxs2 = np.sqrt( (1-fermi(E,mus,betas))*gammas )*ds
    return [auxs1,auxs2]

####################################################
##############aquitmbn##############################
####################################################

def Dd(E,mud,betad,gammad):
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*dddag
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*dd  

    return [auxd1,auxd2]


def Dissipator(E,mus,mud,betas,betad,gammas,gammad):
    DS = Ds(E,mus,betas,gammas)
    DD = Dd(E,mud,betad,gammad)

    tot = []
    for s in DS:
        tot.append(s)
    for d in DD:
        tot.append(d)    

    return tot

def Hamiltonian(E):
    a1 = E*ns + E*nd
    return a1

def Inte(g):
    a2 = g*( np.matmul(dsdag,dd) + np.matmul(dddag,ds) )
    return a2

#########################################################
####################evolutio for time independt things###
#########################################################
#######################(thisisgeneral)###################

def classic(L0,L0draz,P,Q,V,rho0,t):
    #armar PVQLOdrazQVP
    paso1 = np.matmul(V,P)
    paso2 = np.matmul(Q,paso1)
    paso3 = np.matmul(L0draz,paso2)
    paso4 = np.matmul(Q,paso3)
    paso5 = np.matmul(V,paso4)
    paso6 = np.matmul(P,paso5)
    Final = (L0 - paso6)    
    d = len(rho0)
    vec_rho =  np.reshape(rho0,(d**2,1))
    propagator = expm (Final *t)
    prop = P@vec_rho
    vec_rho_t = propagator@prop
    final = P@vec_rho_t

    return np.reshape(final,(d,d))

def ratem(L0,L0draz,P,Q,V):
    #armar PVQLOdrazQVP
    paso1 = np.matmul(V,P)
    paso2 = np.matmul(Q,paso1)
    paso3 = np.matmul(L0draz,paso2)
    paso4 = np.matmul(Q,paso3)
    paso5 = np.matmul(V,paso4)
    paso6 = np.matmul(P,paso5)
    Final = (L0 - paso6)    
     
    return Final 




U0 = 0.
g0 = 0.005

eV = 6.5
mus1 = eV/2
mud1 = -eV/2

betas,betad = 1,1
gs,gd= 1/100,1/100

rho0 = np.array([[1/4,0,0,0],
                 [0,1/4,-0.1j,0],
                 [0,0.1j,1/4,0],
                 [0,0,0,1/4]])

rhovec = Vectorization(rho0)

Ls = Dissipator(0,mus1,mud1,betas,betad,gs,gd)
H = Hamiltonian(0)
L0f = Liouvillian(H,Ls)
#print(len(L0f))
tole = 1E-8
Draz = Drazinalg(L0f,tole)
Draz0 = Drazinspectral(L0f,tole)
P0 = Prin(rho0,basis)
Q0 = Qpart(rho0,P0)


#rhof = Q0@rhovec
#print(vecback(rhof,4))

gss = np.linspace(0.,10,20000)
gaux = []
Il = []
Ilx = []
#I2l = []
#Nss = 3
Wl0 = fermi(0,mus1,betas)*gs
Wdr = fermi(0,mus1,betas)*gs
W0l = (1-fermi(0,mus1,betas))*gs
Wrd = W0l

p00 = []
p10 = []
p01 = []
p11 = []

for g in gss:
    t=10000
    gau = g/gs
    V0 = Inte(g)
    Vnu = pert(V0)
    gaux.append(gau)
    tot = classic(L0f,Draz0,P0,Q0,Vnu,rho0,t)
    W = ratem(L0f,Draz0,P0,Q0,Vnu)
    #print(Il0/gs)  
    #print(g)
    #print(Wl0,W0l)
    #print(W[0,10].real,W[10,0].real)
    p0,pl,pr,pdd = tot[3,3].real,tot[1,1].real,tot[2,2].real,tot[0,0].real
    Il0 = Wl0*p0 - W0l*pl + Wdr*pr - Wrd*pdd
    Ilxx = W[5,15].real*p0 - W[15,5].real*pl + W[0,10].real*pr - W[10,0].real*pdd
    Il.append(Il0/gs)
    Ilx.append(Ilxx/gs)
    #Il.append(tot[0,0].real)
    p00.append(p0)
    p10.append(pl)
    p01.append(pr)
    p11.append(pdd)
    
print(np.kron(v00,v00).T)
print(np.kron(v10,v10).T.real)
print(np.kron(v01,v01).T.real)
print(np.kron(v11,v11).T.real)

plt.plot( gaux,Il)
plt.scatter(gaux,Ilx)
plt.ylabel(r'$I_{L}/\gamma$',fontsize = 20)     
plt.xlabel(r'$g/\gamma$',fontsize = 20)
plt.xscale("log")
plt.show()
#plt.plot( gaux,I2l)
#plt.ylabel(r'$\langle \langle I^{2}_{L} \rangle \rangle/\gamma$',fontsize = 20)     
#plt.xlabel(r'$g/\gamma$',fontsize = 20)
#plt.xscale("log")
#plt.show()


V0f = Inte(0.005)
Vnu0 = pert(V0f)
values =  classic(L0f,Draz,P0,Q0,Vnu0,rho0,8000)

W0 = ratem(L0f,Draz0,P0,Q0,Vnu0)

plt.imshow(values.real)
plt.colorbar()
plt.show()

plt.imshow(W0.real)
plt.colorbar()
plt.show()
    
plt.plot( gaux,p00)
plt.plot( gaux,p10)
plt.plot( gaux,p01)
plt.plot( gaux,p11)
plt.xscale("log")
plt.xlabel(r'$g/\gamma$',fontsize = 20)
plt.show()