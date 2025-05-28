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

nd = np.matmul(dddag,dd)
nl = np.matmul(dldag,dl)
nr = np.matmul(drdag,dr)

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

######################################################################
####################Dissipator3QD#####################################
######################################################################

#operadores del disipador dL
def Dl(E,U,Uf,mul,betal,gammal,gammalU):
    d = len(nl)
    auxl1 = np.sqrt( fermi(E,mul,betal)*gammal )*np.matmul( np.matmul((np.eye(d)-nr),(np.eye(d)-nd)),dldag )
    auxl2 = np.sqrt( (1-fermi(E,mul,betal))*gammal )*np.matmul( np.matmul((np.eye(d)-nr),(np.eye(d)-nd)),dl)
    auxl3 = np.sqrt( fermi(E+U,mul,betal)*gammalU )*np.matmul( np.matmul((np.eye(d)-nr) ,nd),dldag )
    auxl4 = np.sqrt( (1-fermi(E+U,mul,betal))*gammalU )*np.matmul(np.matmul((np.eye(d)-nr) ,nd),dl)
    auxl5 = np.sqrt( fermi(E+Uf,mul,betal)*gammal )*np.matmul( np.matmul((np.eye(d)-nd) ,nr),dldag )
    auxl6 = np.sqrt( (1-fermi(E+Uf,mul,betal))*gammal )*np.matmul(np.matmul((np.eye(d)-nd) ,nr),dl)
    auxl7 = np.sqrt( fermi(E+U+Uf,mul,betal)*gammal )*np.matmul( np.matmul(nr,nd),dldag )
    auxl8 = np.sqrt( (1-fermi(E+U+Uf,mul,betal))*gammal )*np.matmul(np.matmul(nr,nd),dl)

    return [auxl1,auxl2,auxl3,auxl4,auxl5,auxl6,auxl7,auxl8]
#operadores del disipador dr
def Dr(E,U,Uf,mur,betar,gammar,gammarU):
    d = len(nr)
    auxr1 = np.sqrt( fermi(E,mur,betar)*gammar )*np.matmul( np.matmul((np.eye(d)-nl),(np.eye(d)-nd)),drdag )
    auxr2 = np.sqrt( (1-fermi(E,mur,betar))*gammar )*np.matmul( np.matmul((np.eye(d)-nl),(np.eye(d)-nd)),dr)
    auxr3 = np.sqrt( fermi(E+U,mur,betar)*gammarU )*np.matmul( np.matmul((np.eye(d)-nl) ,nd),drdag )
    auxr4 = np.sqrt( (1-fermi(E+U,mur,betar))*gammarU )*np.matmul(np.matmul((np.eye(d)-nl) ,nd),dr)
    auxr5 = np.sqrt( fermi(E+Uf,mur,betar)*gammar )*np.matmul( np.matmul((np.eye(d)-nd) ,nl),drdag )
    auxr6 = np.sqrt( (1-fermi(E+Uf,mur,betar))*gammar )*np.matmul(np.matmul((np.eye(d)-nd) ,nl),dr)
    auxr7 = np.sqrt( fermi(E+U+Uf,mur,betar)*gammar )*np.matmul( np.matmul( nl,nd),drdag )
    auxr8 = np.sqrt( (1-fermi(E+U+Uf,mur,betar))*gammar )*np.matmul(np.matmul(nl ,nd),dr)

    return [auxr1,auxr2,auxr3,auxr4,auxr5,auxr6,auxr7,auxr8]

#operadores del disipador dd
def Dd(E,U,mud,betad,gammad,gammadU):
    d = len(nr)
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*np.matmul( np.matmul((np.eye(d)-nl),(np.eye(d)-nr)),dddag )
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*np.matmul( np.matmul((np.eye(d)-nl),(np.eye(d)-nr)),dd)
    auxd3 = np.sqrt( fermi(E+U,mud,betad)*gammadU )*np.matmul( np.matmul((np.eye(d)-nl) ,nr) + np.matmul((np.eye(d)-nr) ,nl) ,dddag )
    auxd4 = np.sqrt( (1-fermi(E+U,mud,betad))*gammadU )*np.matmul(np.matmul((np.eye(d)-nl) ,nr) + np.matmul((np.eye(d)-nr) ,nl),dd)
    auxd5 = np.sqrt( fermi(E+ (2*U),mud,betad)*gammadU )*np.matmul( np.matmul(nr ,nl),dddag )
    auxd6 = np.sqrt( (1-fermi(E+(2*U),mud,betad))*gammadU )*np.matmul(np.matmul(nr ,nl),dd)

    return [auxd1,auxd2,auxd3,auxd4,auxd5,auxd6]

def Dissipator(E,Ed,U,Uf,mul,mur,mud,betal,betar,betad,gammal,gammalU,gammar,gammarU,gammad,gammadU):
    DR = Dr(E,U,Uf,mur,betar,gammar,gammarU)
    DL = Dl(E,U,Uf,mul,betal,gammal,gammalU)
    DD = Dd(Ed,U,mud,betad,gammad,gammadU)

    tot = []
    for l in DL:
        tot.append(l)
    for r in DR:
        tot.append(r)
    for d in DD:
        tot.append(d)    

    return tot

def Htd(E,Ed,U,Uf):
    Htdl = E*nl + E*nr + Ed*nd + U*(np.matmul(nr,nd) + np.matmul(nl,nd)) + Uf*np.matmul(nr,nl) 
    return Htdl

def Inte(g):
    a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
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

#######################################
##############base#####################
#p000,p100,p010,p001,p110,p011,p101,p111
##################################
############ratematrix############
##################################
def ratel(W0l,Wl0,Wluf,Wufl,Wlu,Wul,Wu2l,Wlu2):
    ####base###

    Wl = np.zeros((8,8))
    #W0l
    Wl[0,0],Wl[0,1],Wl[1,0],Wl[1,1] = -Wl0,W0l,Wl0,-W0l
    #Wufl
    Wl[2,2],Wl[2,4],Wl[4,2],Wl[4,4] = -Wufl,Wluf,Wufl,-Wluf
    #Wul
    Wl[3,3],Wl[3,6],Wl[6,3],Wl[6,6] = -Wul,Wlu,Wul,-Wlu
    #Wu2l
    Wl[5,5],Wl[5,7],Wl[7,5],Wl[7,7] = -Wu2l,Wlu2,Wu2l,-Wlu2
    return Wl

def rater(W0r,Wr0,Wruf,Wufr,Wru,Wur,Wu2r,Wru2):
    ####base###

    Wr = np.zeros((8,8))
    #W0r
    Wr[0,0],Wr[0,2],Wr[2,0],Wr[2,2] = -Wr0,W0r,Wr0,-W0r
    #Wufr
    Wr[1,1],Wr[1,4],Wr[4,1],Wr[4,4] = -Wufr,Wruf,Wufr,-Wruf
    #Wur
    Wr[3,3],Wr[3,5],Wr[5,3],Wr[5,5] = -Wur,Wru,Wur,-Wru
    #Wu2r
    Wr[6,6],Wr[6,7],Wr[7,6],Wr[7,7] = -Wu2r,Wru2,Wu2r,-Wru2
    return Wr

def rateg(Wlr,Wrl,Wlru,Wrlu):
    Wcoup = np.zeros((8,8))
    #Wlr
    Wcoup[1,1],Wcoup[1,2],Wcoup[2,1],Wcoup[2,2] = -Wrl,Wlr,Wrl,-Wlr
    #Wlru
    Wcoup[7,7],Wcoup[7,6],Wcoup[6,7],Wcoup[6,6] = -Wrlu,Wlru,Wrlu,-Wlru
    return Wcoup

def ratelr(Wl,Wr,Wlr):

    return Wl+Wr+Wlr

def vecflow(W,p,Jn):
    return Jn.T@W@p



betar,betad,betal = 1/100,1/2,1/100

gr,grU = (1/100)*(1/6), 1/100
gl,glU = 1/100, (1/100)*(1/6)
gd,gdU = 1/50,1/50


rho0 = np.array([[1/8,0,0,0,0,0,0,0],
                 [0,1/8,0,0,0,0,0,0],
                 [0,0,1/8,0,0,0,0,0],
                 [0,0,0,1/8,0,0,0,0],
                 [0,0,0,0,1/8,0,0,0],
                 [0,0,0,0,0,1/8,0,0],
                 [0,0,0,0,0,0,1/8,0],
                 [0,0,0,0,0,0,0,1/8]])

Num = 200
eVs0 = np.linspace(0,800,Num)
g0 = 5/1000

print("p000",v00.T)
print("p100",v100.T.real)
print("p010",v010.T.real)
print("p001",v001.T.real)
print("p110",v110.T.real)
print("p011",v011.T.real)
print("p101",v101.T.real)
print("p111",v111.T.real)

If = []
eVs = []
curr = []
Probnt1f = []
Probnt2f = []
Probnt3f = []
Probnt4f = []
Probnt5f = []
Probnt6f = []
Probnt7f = []
Probnt8f = []
for ev in eVs0:
    mud0 = 2
    U00 = 40 #10
    #mud0 = 1-U00/2
    #Con Ed0 = mud0 -U00/2,E0=4 hay flujo de energia pero un orden menor al de
    #flujo de informacion
    #Ed0 = 1
    Ed0 = mud0 -U00/2
    Uf0 = 500 #50
    #Probar condicion (U00/E0)<<1,Strasberg
    E0 = 0
    Ls0 = Dissipator(E0,Ed0,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,gl,glU,gr,grU,gd,gdU)
    H0 = Htd(E0,Ed0,U00,Uf0)
    L0 = Liouvillian(H0,Ls0)
    t=10000    
    tole = 1E-6
    Draz0 = Drazinspectral(L0,tole)
    P0 = Prin(rho0,basis)
    Q0 = Qpart(rho0,P0)
    V0 = Inte(g0)
    Vnu = pert(V0)
    tot = classic(L0,Draz0,P0,Q0,Vnu,rho0,t)
    W = ratem(L0,Draz0,P0,Q0,Vnu)
    p000,p100,p010,p001 = tot[7,7].real,tot[3,3].real,tot[5,5].real,tot[6,6].real
    p110,p011,p101 = tot[1,1].real,tot[4,4].real,tot[2,2].real
    p111 = tot[0,0].real
    ps = np.array([[p000],
                [p100],
                [p010],
                [p001],
                [p110],
                [p011],
                [p101],
                [p111]])
    #informationvecto
    logps = np.array([[np.log(p000)],
                [np.log(p100)],
                [np.log(p010)],
                [np.log(p001)],
                [np.log(p110)],
                [np.log(p011)],
                [np.log(p101)],
                [np.log(p111)]])
    
    #numbervector 
    N = np.array([[0],
                [1],
                [1],
                [1],
                [2],
                [2],
                [2],
                [3]])
    
    W0l,Wl0,Wluf,Wufl = W[63,27].real,W[27,63].real,W[45,9].real,W[9,45].real
    Wlu,Wul,Wlu2,Wu2l = W[54,18].real,W[18,54].real,W[36,0].real,W[0,36].real

    W0r,Wr0,Wruf,Wufr = W[63,45].real,W[45,63].real,W[27,9].real,W[9,27].real
    Wru,Wur,Wru2,Wu2r = W[54,36].real,W[36,54].real,W[18,0].real,W[0,18].real

    Wlr,Wrl,Wlru,Wrlu = W[27,45].real,W[45,27].real,W[18,36].real,W[36,18].real
    Probnt1f.append(p111)
    Probnt2f.append(p110)
    Probnt3f.append(p101)
    Probnt4f.append(p100)
    Probnt5f.append(p011)
    Probnt6f.append(p010)
    Probnt7f.append(p001)
    Probnt8f.append(p000)
    Wl = ratel(W0l,Wl0,Wluf,Wufl,Wlu,Wul,Wu2l,Wlu2)
    Wr = rater(W0r,Wr0,Wruf,Wufr,Wru,Wur,Wu2r,Wru2)
    WLR = rateg(Wlr,Wrl,Wlru,Wrlu)
    Wfs = Wl+WLR
    Wf = ratelr(Wl,Wr,WLR)

    curre = vecflow(Wfs,ps,N)
    infos = vecflow(Wf,ps,logps)
    eVs.append(ev*betal)
    If.append(infos[0][0])
    curr.append(curre[0][0])



plt.plot(eVs,Probnt1f,linestyle='--', dashes=(5, 9), color='green',lw = 4,label = r'$\rho_{111}$')
plt.plot(eVs,Probnt2f, color='green',lw = 4,label = r'$\rho_{110}$')
plt.plot(eVs,Probnt3f,linestyle='--', dashes=(5, 9), color='black',lw=3,label = r'$\rho_{101}$')
plt.plot(eVs,Probnt4f, color='black',lw = 4, label = r'$\rho_{100}$')
plt.plot(eVs,Probnt5f,color = 'orange',lw = 4, label = r'$\rho_{011}$')
plt.plot(eVs,Probnt6f, color='red',lw = 4, label = r'$\rho_{010}$')
plt.plot(eVs,Probnt7f, color='blue',lw = 4, label = r'$\rho_{001}$')
plt.plot(eVs,Probnt8f,color = 'm', lw = 4,label = r'$\rho_{000}$')
plt.legend(loc = "upper right",fontsize=15)
plt.xlabel(r'$eV/T$',fontsize=25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.show()

plt.plot(eVs,curr,label = r'$\dot{N}_{L}$', color = 'b')
#plt.xscale("log")
plt.xlabel(r'$eV$',fontsize = 20)
#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  
#plt.ylabel("Particle current",fontsize = 20)
plt.legend()
plt.show()             

plt.plot(eVs,If,label = r'$\dot{I}_{LR}$', color = 'b')
#plt.xscale("log")
plt.xlabel(r'$eV$',fontsize = 20)
#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  
#plt.ylabel("Particle current",fontsize = 20)
plt.legend()
plt.xscale("log")
plt.show()   


archivo = open("classic","w")
decimal_places = 7
total_width = 8
format_str = f"{{:.{decimal_places}f}}" 
#format_str = f"{{:{total_width}.{decimal_places}f}}"
for i in range(Num):
    #print(curr[i])
    archivo.write( format_str.format(eVs[i])) #guarda el grado del nodo
    #archivo.write(str(xs[i])) 
    archivo.write(" ") 
    #archivo.write(str(ys[i]))
    archivo.write( format_str.format(-curr[i]))
    archivo.write(" ")  
    archivo.write( format_str.format(Probnt1f[i]))
    archivo.write(" ")
    archivo.write( format_str.format(Probnt2f[i]))
    archivo.write(" ")
    archivo.write( format_str.format(Probnt3f[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(Probnt4f[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(Probnt5f[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(Probnt6f[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(Probnt7f[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(Probnt8f[i]))
    archivo.write(" ") 
    #archivo.write(str(ys[i]))
    archivo.write( format_str.format(-If[i]))
    archivo.write("\n")

