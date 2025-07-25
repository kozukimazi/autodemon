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

def conmutator(A,B):
    return np.matmul(A,B) - np.matmul(B,A)
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

#this is what we use
def partial_trace_d(A):
    """
    Partial trace over the third qubit (subsystem d).
    Input:  A — 8x8 matrix (operator on 3 qubits)
    Output: 4x4 matrix (operator on qubits l and r)
    """
    #A = A_{ijk;i'j'k'}|ijk><i'j'k'|
    #(this is default numpy)
    #the reshape is for rows i = 4l+2r+d
    #for columns is j = 4l' + 2r' + d'
    A = A.reshape(2, 2, 2, 2, 2, 2)  # reshape to (l, r, d, l', r', d')
    return np.trace(A, axis1=2, axis2=5).reshape(4,4)  # trace over d and d'

# d is MSB, then r, then l
def partial_trace_d_custom(A):
    A = A.reshape(2, 2, 2, 2, 2, 2)  # (d, r, l, d', r', l')
    # trace out d and d' → axes 0 and 3
    return np.trace(A, axis1=0, axis2=3).reshape(4,4)


def Liouvillian( H,Ls, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    superL = sum( [np.kron(L.conjugate(),L) - 1/2 * (np.kron( np.eye(d), L.conjugate().T.dot(L)) +
                                                     np.kron( L.T.dot(L.conjugate()),np.eye(d) ))
                                                      for L in Ls ] )        
    return superH + superL

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

#Funcion auxiliar para calcular los D_{i}\rho

def Dissipate(H,Ls, rho):
    d = len(H)
    superL = np.zeros((d,d) , dtype = np.complex_)
    for L in Ls:
        d = np.matmul( np.matmul( L,rho ),L.transpose().conjugate() )
        e = (1/2)*anticonmutador( np.matmul(L.transpose().conjugate(),L), rho )
        superL += d-e
        
    return superL

def Liouvillian( H,Ls, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    superL = sum( [np.kron(L.conjugate(),L) - 1/2 * (np.kron( np.eye(d), L.conjugate().T.dot(L)) +
                                                     np.kron( L.T.dot(L.conjugate()),np.eye(d) ))
                                                      for L in Ls ] )    
    return superH + superL

def Propagate(rho0,superop,t):
    d = len(rho0)
    propagator = expm (superop *t)
    vec_rho_t = propagator @ np.reshape(rho0,(d**2,1))
    return np.reshape( vec_rho_t, (d,d) )

def Hamiltonian(E,Ed,U,Uf,g):
    a1 = E*nl + E*nr + Ed*nd
    a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(nl,nd) +  np.matmul(nr,nd) ) + Uf*np.matmul(nl,nr) 
    return a1+a2+a3

def Htd(E,Ed,U,Uf):
    Htdl = E*nl + E*nr + Ed*nd + U*(np.matmul(nr,nd) + np.matmul(nl,nd)) + Uf*np.matmul(nr,nl) 
    return Htdl

def logaux(x):
    a= 0
    if (abs(x)> 0.000001):
       a += np.log(x)

    return a    

def arctanaux(x,y):

    if (abs(y)<1E-7):
        return np.pi/2
    
    else:
        return np.arctan(x/y)


def currentsnew(Htd,Hs,mul,mur,mud,Ll,Lr,Ld,superop,rho0,t):
    Nop = nl + nr + nd
    rhof = Propagate(rho0,superop,t)
    cohe = abs(rhof[5,3]) + abs(rhof[4,2]) 
    cohesum = abs(rhof[5,3] + rhof[4,2])
    PD = rhof[0,0].real + rhof[1,1].real 
    P0 = rhof[7,7].real + rhof[6,6].real 
    concurrencef = 2*cohesum - 2*np.sqrt(P0*PD)
    concf = 0
    if (concurrencef > 0):
        concf += concurrencef
    else: 
        concf = 0
    Dl = Dissipate(Htd,Ll,rhof)
    Dr = Dissipate(Htd,Lr,rhof)
    Dd = Dissipate(Htd,Ld,rhof)

    Dt = Dl + Dr + Dd
    rhofLR = partial_trace_d(rhof)
    rhofLRU = -partial_trace_d(np.matmul(Dt,rhof) )
    #rhofLRU = -1j*partial_trace_d( conmutator(Hs,rhof) )
    rhofLRD = partial_trace_d(np.matmul(Dt,rhof) )

    Qopl = (Htd - mul*Nop)
    Qopr = (Htd - mur*Nop)
    Qopd = (Htd - mud*Nop)

    aux = logm(rhof)
    auxlr = logm(rhofLR)
    Ql = np.trace( np.matmul( Dl,Qopl  ) )
    Qr = np.trace( np.matmul( Dr,Qopr  ) )
    Qd = np.trace( np.matmul( Dd,Qopd  ) )

    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )

    IU = -np.trace( np.matmul( rhofLRU,auxlr ) )
    Idiss = -np.trace( np.matmul( rhofLRD,auxlr ) ) -Sl-Sr

    El = np.trace( np.matmul( Dl,Htd  ) )
    Er = np.trace( np.matmul( Dr,Htd  ) )
    Ed = np.trace( np.matmul(Dd,Htd))

    Wl = mul*np.trace(np.matmul( Dl, Nop ))
    Wr = mur*np.trace(np.matmul( Dr, Nop ))
    Wd = mud*np.trace(np.matmul( Dd, Nop ))

    Nl = np.trace( np.matmul(Dl,Nop ) )
    Nr = np.trace(np.matmul(Dr,Nop ))
    Nd = np.trace(np.matmul( Dd,Nop ))
    return Ql.real, Qr.real, Qd.real, Sl.real, Sr.real,Sd.real, El.real, Er.real, Ed.real, Wl.real,Wr.real,Wd.real,cohe,concf,Nl.real,Nr.real,Nd.real,IU.real,Idiss.real


E = 0
U0 = 40.
Uf = 500
g0 = 5/1000
#g0 = 0.000001
eV = 450
mul1 = eV/2
mur1 = -eV/2
mud1 = 2
Ed = mud1-U0/2
#betar,betad,betal = 1/50,1/20,1/50
#comportamiento raro
#betar,betad,betal = 1/10,1/5,1/10
#demonio
#betar,betad,betal = 1/100,5/10,1/100
#refrigerador 
#betar,betad,betal = 1/10,20,1/10
#strasberg
#betar,betad,betal = 1/100,10,1/100
#transport
betar,betad,betal = 1/100,1/2,1/100

n=6
gr,grU = (1/100)*(1/n), 1/100
gl,glU = 1/100, (1/100)*(1/n)
gd,gdU = 1/50,1/50 


alp,alp2,alp3,alp4,alp5 = 0.,0.,0.0,0.0,0.
a,b,c,d = 1j*alp,1j*alp2,1j*alp3,1j*alp4
#111,110,101,100,011,010,001,000
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
Ql = []
Qr = []
Qd = []
Sls = []
Srs = []
Sds = []
Slr = []
entropf = []
Isl = []
Id = []
Els = []
Ers = []
Eds = []
Erl = []
Qlr = []
Wl = []
Wr = []
Wd = []
Wt = []
Flr = []
Tisl = []
cohes = []
Wdf = []
Qdf = []
Tid = []
Fd = []
Ilf = []
Irf = []
Nls = []
Nrs = []
Nds = []
cohev = []
concv = []
eVs = []
auxff = []
Probnt10 = []
Probnt20 = []
Probnt30 = []
Probnt40 = []
Probnt50 = []
Probnt60 = []
Probnt70 = []
Probnt80 = []
resta = []
cohesex = []
infolrU = []
infolrdis = []
infosum = []
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
    H0 = Hamiltonian(E0,Ed0,U00,Uf0,g0)
    superop0 = Liouvillian(H0,Ls0)
    Ll0 = Dl(E0,U00,Uf0,ev/2,betal,gl,glU)
    Lr0 = Dr(E0,U00,Uf0,-ev/2,betar,gr,grU)
    Ld0 = Dd(Ed0,U00,mud0,betad,gd,gdU)
    Htd0 =  Htd(E0,Ed0,U00,Uf0)
    t=40000
    rhof = Propagate(rho0,superop0,t)
    Ql0,Qr0,Qd0,Sl0,Sr0,Sd0,El0,Er0,Ed0,Wl0,Wr0,Wd0,cohe0,concu0,Nl0,Nr0,Nd0,ilru,ilrdis = currentsnew(Htd0,H0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,superop0,rho0,t)
    Ql.append(Ql0)
    Qr.append(Qr0)
    Qd.append(Qd0)
    #cohev.append(abs(rhof[5,3]) + abs(rhof[4,2]) )
    concv.append(concu0)    
    sigmal = Sl0 - betal*Ql0
    sigmar = Sr0 - betar*Qr0
    Sls.append(Sl0 - betal*Ql0)
    Srs.append(Sr0 - betar*Qr0)
    Sds.append(Sd0 - betad*Qd0)
    Slr.append( sigmal + sigmar )
    Isl0 = -Sl0 - Sr0
    Id0 = -Sd0
    Il0 = -Sl0
    Ir0 = -Sr0
    Isl.append(Isl0)
    Id.append(-Sd0)
    Els.append(El0)
    Ers.append(Er0)
    Eds.append(Ed0)
    Erl.append(El0 + Er0 )
    Qlr.append(Ql0+Qr0 )
    Flr.append(El0 + Er0 + (1/betal)*Isl0 )
    Tisl.append((1/betal)*Isl0  )
    Wt.append(Wl0 + Wr0 )
    Wdf.append(Wd0)
    Qdf.append(Qd0)
    Fd.append(Ed0 + (1/betad)*Id0 )
    Tid.append((1/betad)*Id0  )
    cohes.append(cohe0)
    Ilf.append(-Sl0)
    Irf.append(-Sr0)
    Nls.append(Nl0)
    Nrs.append(Nr0)
    Nds.append(Nd0)
    eVs.append(ev*betal)
    entropf.append( -betal*(Ql0+Qr0) )
    auxff.append(0)
    Probnt10.append(rhof[0,0].real )
    Probnt20.append(rhof[1,1].real )
    Probnt30.append(rhof[2,2].real )
    Probnt40.append(rhof[3,3].real )
    Probnt50.append(rhof[4,4].real )
    Probnt60.append(rhof[5,5].real )
    Probnt70.append(rhof[6,6].real )
    Probnt80.append(rhof[7,7].real )
    infolrU.append(ilru.real)
    infolrdis.append(ilrdis.real)
    infosum.append(ilrdis.real + ilru.real)
    resta.append(abs( rhof[3,3].real - rhof[5,5].real ))
    cohesex.append(abs( rhof[3,5]).real)


plt.plot(eVs,infolrU,linestyle='--', dashes=(5, 9), color='red',lw = 4,label = r'$\dot{I}^{U}_{LR}$')
plt.plot(eVs,infolrdis,linestyle='--', dashes=(5, 9), color='blue', lw=4,label = r'$\dot{I}^{D}_{LR}$') 
#plt.plot(eVs,infosum,linestyle='--', dashes=(5, 9), color='black',lw=4,label = r'$\dot{I}_{LR}$')
#plt.plot(eVs,Isl,linestyle='--', dashes=(5, 9), color='green',lw=4,label = r'$\dot{I}_{LRori}$')
#plt.plot(eVs,Nls,label = r'$\dot{N}_{L}$')
#plt.plot(eVs,Nrs, label = r'$\dot{N}_{R}$') 
#plt.plot(eVs,Nds,label = r'$\dot{N}_{d}$')
#plt.xscale("log")
plt.xlabel(r'$eV/T$',fontsize = 20)
#plt.ylabel(r'$J_{\alpha}$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)  # Y-axis tick labels
plt.legend(loc = "lower left",fontsize=15) 
plt.show()   



