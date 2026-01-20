import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import cmath


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


nd = np.matmul(dddag,dd)

##definir base con objetivo de calcular coherencias
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


def twolevel(el,er,g):
    Delta = (el-er)/2
    if (Delta**2 + g**2>1E-9):
        aux = Delta/np.sqrt( Delta**2 + g**2)
        theta = math.acos(aux)
    #revisar este else
    else:
        theta = math.acos(1)    
    dup = np.cos(theta/2)*dl + np.sin(theta/2)*dr
    dmin =  np.cos(theta/2)*dr - np.sin(theta/2)*dl
    dupdag = np.cos(theta/2)*dldag + np.sin(theta/2)*drdag
    dmindag =  np.cos(theta/2)*drdag - np.sin(theta/2)*dldag 
    return dup,dupdag,dmin,dmindag

def energy(el,er,g):
    Delta = (el-er)/2
    if (Delta**2 + g**2>1E-9):
        aux = Delta/np.sqrt( Delta**2 + g**2)
        theta = math.acos(aux)
    else:
        theta = math.acos(1) 

    return Delta, g, theta


def Liouvillian( H,Ls, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    superL = sum( [np.kron(L.conjugate(),L) - 1/2 * (np.kron( np.eye(d), L.conjugate().T.dot(L)) +
                                                     np.kron( L.T.dot(L.conjugate()),np.eye(d) ))
                                                      for L in Ls ] )        
    return superH + superL

#operadores del disipador dd
def Dd(El,Er,g,E,U,mud,betad,gammad):
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    nm = np.matmul(dmindag,dmin)
    np0 = np.matmul(dupdag,dup)
    d = 8
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nm)),dddag )
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nm)),dd)
    auxd3 = np.sqrt( fermi(E+U,mud,betad)*gammad )*np.matmul( np.matmul((np.eye(d)-np0) ,nm) + np.matmul((np.eye(d)-nm) ,np0) ,dddag )
    auxd4 = np.sqrt( (1-fermi(E+U,mud,betad))*gammad )*np.matmul(np.matmul((np.eye(d)-np0) ,nm) + np.matmul((np.eye(d)-nm) ,np0),dd)
    auxd5 = np.sqrt( fermi(E+ (2*U),mud,betad)*gammad )*np.matmul( np.matmul(np0 ,nm),dddag )
    auxd6 = np.sqrt( (1-fermi(E+(2*U),mud,betad))*gammad )*np.matmul(np.matmul(np0 ,nm),dd)

    return [auxd1,auxd2,auxd3,auxd4,auxd5,auxd6]

#operadores del disipador positivos
def Dp(Eup,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dup,dupdag,dmin,dmindag,theta):
    d = 8
    Ep = Eup
    nm = np.matmul(dmindag,dmin)
    Cl = np.cos(theta/2)
    Cr = np.sin(theta/2)
    auxl1 = np.sqrt( fermi(Ep,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dupdag )
    auxl2 = np.sqrt( (1-fermi(Ep,mul,betal))*gl )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dup)
    auxl3 = np.sqrt( fermi(Ep+U,mul,betal)*glU )*np.matmul( np.matmul((np.eye(d)-nm) ,nd),dupdag )
    auxl4 = np.sqrt( (1-fermi(Ep+U,mul,betal))*glU )*np.matmul(np.matmul((np.eye(d)-nm) ,nd),dup)
    auxl5 = np.sqrt( fermi(Ep+Uf,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nd) ,nm),dupdag )
    auxl6 = np.sqrt( (1-fermi(Ep+Uf,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-nd) ,nm),dup)
    auxl7 = np.sqrt( fermi(Ep+U+Uf,mul,betal)*gl )*np.matmul( np.matmul(nm,nd),dupdag )
    auxl8 = np.sqrt( (1-fermi(Ep+U+Uf,mul,betal))*gl )*np.matmul(np.matmul(nm,nd),dup)

    auxr1 = np.sqrt( fermi(Ep,mur,betar)*gr )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dupdag )
    auxr2 = np.sqrt( (1-fermi(Ep,mur,betar))*gr )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dup)
    auxr3 = np.sqrt( fermi(Ep+U,mur,betar)*grU )*np.matmul( np.matmul((np.eye(d)-nm) ,nd),dupdag )
    auxr4 = np.sqrt( (1-fermi(Ep+U,mur,betar))*grU )*np.matmul(np.matmul((np.eye(d)-nm) ,nd),dup)
    auxr5 = np.sqrt( fermi(Ep+Uf,mur,betar)*gr )*np.matmul( np.matmul((np.eye(d)-nd) ,nm),dupdag )
    auxr6 = np.sqrt( (1-fermi(Ep+Uf,mur,betar))*gr )*np.matmul(np.matmul((np.eye(d)-nd) ,nm),dup)
    auxr7 = np.sqrt( fermi(Ep+U+Uf,mur,betar)*gr )*np.matmul( np.matmul(nm,nd),dupdag )
    auxr8 = np.sqrt( (1-fermi(Ep+U+Uf,mur,betar))*gr )*np.matmul(np.matmul(nm,nd),dup)

        

    Dpl = [Cl*auxl1, Cl*auxl2 , Cl*auxl3, Cl*auxl4, Cl*auxl5, Cl*auxl6, Cl*auxl7, Cl*auxl8]
    Dpr = [Cr*auxr1, Cr*auxr2 , Cr*auxr3, Cr*auxr4, Cr*auxr5, Cr*auxr6, Cr*auxr7, Cr*auxr8]
    
    return Dpl,Dpr

#operadores del disipador dpositivos
def Dm(Emin,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dmin,dmindag,dup,dupdag,theta):
    d = 8
    Em = Emin
    np0 = np.matmul(dupdag,dup)
    Cl = -np.sin(theta/2)
    Cr = np.cos(theta/2)
    auxl1 = np.sqrt( fermi(Em,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dmindag )
    auxl2 = np.sqrt( (1-fermi(Em,mul,betal))*gl )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dmin)
    auxl3 = np.sqrt( fermi(Em+U,mul,betal)*glU )*np.matmul( np.matmul((np.eye(d)-np0) ,nd),dmindag )
    auxl4 = np.sqrt( (1-fermi(Em+U,mul,betal))*glU )*np.matmul(np.matmul((np.eye(d)-np0) ,nd),dmin)
    auxl5 = np.sqrt( fermi(Em+Uf,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nd) ,np0),dmindag )
    auxl6 = np.sqrt( (1-fermi(Em+Uf,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-nd) ,np0),dmin)
    auxl7 = np.sqrt( fermi(Em+U+Uf,mul,betal)*gl )*np.matmul( np.matmul(np0,nd),dmindag )
    auxl8 = np.sqrt( (1-fermi(Em+U+Uf,mul,betal))*gl )*np.matmul(np.matmul(np0,nd),dmin)

    auxr1 = np.sqrt( fermi(Em,mur,betar)*gr )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dmindag )
    auxr2 = np.sqrt( (1-fermi(Em,mur,betar))*gr )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dmin)
    auxr3 = np.sqrt( fermi(Em+U,mur,betar)*grU )*np.matmul( np.matmul((np.eye(d)-np0) ,nd),dmindag )
    auxr4 = np.sqrt( (1-fermi(Em+U,mur,betar))*grU )*np.matmul(np.matmul((np.eye(d)-np0) ,nd),dmin)
    auxr5 = np.sqrt( fermi(Em+Uf,mur,betar)*gr )*np.matmul( np.matmul((np.eye(d)-nd) ,np0),dmindag )
    auxr6 = np.sqrt( (1-fermi(Em+Uf,mur,betar))*gr )*np.matmul(np.matmul((np.eye(d)-nd) ,np0),dmin)
    auxr7 = np.sqrt( fermi(Em+U+Uf,mur,betar)*gr )*np.matmul( np.matmul(np0,nd),dmindag )
    auxr8 = np.sqrt( (1-fermi(Em+U+Uf,mur,betar))*gr )*np.matmul(np.matmul(np0,nd),dmin)

    Dml = [Cl*auxl1, Cl*auxl2 , Cl*auxl3, Cl*auxl4, Cl*auxl5, Cl*auxl6, Cl*auxl7, Cl*auxl8]
    Dmr = [Cr*auxr1, Cr*auxr2 , Cr*auxr3, Cr*auxr4, Cr*auxr5, Cr*auxr6, Cr*auxr7, Cr*auxr8]
    
    return Dml, Dmr

#####################################################
############phononsection############################
#####################################################
def bose(omega,betaph):

    n = 1/(np.exp(omega*betaph) - 1)
    return n

def Jph(omega,omegac,J0):
    if (omega < 0):
        return -J0*omega*np.exp(-abs(omega)/omegac)
    else:
        return J0*omega*np.exp(-abs(omega)/omegac)

def gammaph(omega,omegac,J0,betaph):
    #omegac = 1E-2
    #gamma0 = 1
    #J0 = 2
    
    if (abs(omega) > 1E-7):
        
        n = bose(omega,betaph)
        #print("aqui tmbn")
        if (omega > 0):
        
            rate =  Jph(omega,omegac,J0)*(n+1)
        else:
            rate = abs(Jph(abs(omega),omegac,J0)*n)
        return rate
        
    else:
        #print("valor:")
        return J0/betaph 

def Dphonon(El,Er,g,betaph,J0,omegac):
    Delta, g, theta = energy(El,Er,g)
    eps = (El+Er)/2
    Eup, Emin = eps + np.sqrt(Delta**2 + g**2), eps - np.sqrt(Delta**2 + g**2)
    
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    nup = np.matmul(dupdag,dup)
    nmin = np.matmul(dmindag,dmin)
    #print("la energia E- - E+ es:")
    #print(Emin-Eup)
    #print(gammaph(Emin-Eup,omegac,J0,betaph))
    auxl1 = np.cos(theta)*np.sqrt( gammaph(Eup-Emin,omegac,J0,betaph) )*np.matmul( dupdag,dmin )
    auxl2 = np.cos(theta)*np.sqrt( gammaph(Emin-Eup,omegac,J0,betaph) )*np.matmul( dmindag,dup )
    auxl3 = np.sin(theta)*np.sqrt( gammaph(0,omegac,J0,betaph) )*(nup-nmin) 
    return [auxl1,auxl2,auxl3]


def Dissipator(El,Er,g,Ed,U,Uf,mul,mur,mud,betal,betar,betad,betaph,gl,glU,gr,grU,gammad,J0,omegac):
    Delta, g, theta = energy(El,Er,g)
    eps = (El+Er)/2
    Eup, Emin = eps + np.sqrt(Delta**2 + g**2), eps - np.sqrt(Delta**2 + g**2)
    
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    DLP,DRP = Dp(Eup,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dup,dupdag,dmin,dmindag,theta)
    DLM,DRM = Dm(Emin,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dmin,dmindag,dup,dupdag,theta)
   
    DD = Dd(El,Er,g,Ed,U,mud,betad,gammad)
    Dph = Dphonon(El,Er,g,betaph,J0,omegac)

    tot = []
    for l in DLP:
        tot.append(l)
    for l in DLM:
        tot.append(l)
    for r in DRP:
        tot.append(r)
    for r in DRM:
        tot.append(r)
    for d in DD:
        tot.append(d)    
    for dph in Dph:
        tot.append(dph)    

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

def Hamiltonian(El,Er,Ed,U,Uf,g):
    Delta, g, theta = energy(El,Er,g)
    eps = (El+Er)/2
    Eup, Emin = eps + np.sqrt(Delta**2 + g**2), eps - np.sqrt(Delta**2 + g**2)
    
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    nm = np.matmul(dmindag,dmin)
    np0 = np.matmul(dupdag,dup)
    #lo diagonalizo aqui
    a1 = Eup*np0 + Emin*nm + Ed*nd
    #a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(np0,nd) +  np.matmul(nm,nd) ) + Uf*np.matmul(np0,nm) 
    return a1+a3

def currents(El,Er,g,Hs,mul,mur,mud,Ll,Lr,Ld,Dph,superop,rho0,t):
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    nm = np.matmul(dmindag,dmin)
    np0 = np.matmul(dupdag,dup)
    Nop = np0 + nm + nd
    rhof = Propagate(rho0,superop,t)

    Dl = Dissipate(Hs,Ll,rhof)
    Dr = Dissipate(Hs,Lr,rhof)
    Dd = Dissipate(Hs,Ld,rhof)
    Dphf = Dissipate(Hs,Dph,rhof)
    Qopl = (Hs - mul*Nop)
    Qopr = (Hs - mur*Nop)
    Qopd = (Hs - mud*Nop)

    
    aux = logm(rhof)
    Ql = np.trace( np.matmul( Dl,Qopl  ) )
    Qr = np.trace( np.matmul( Dr,Qopr  ) )
    Qd = np.trace( np.matmul( Dd,Qopd  ) )
    Qphlr = np.trace( np.matmul( Dphf,Hs  ) )


    Nl = np.trace( np.matmul( Dl,Nop  ) )
    Nr = np.trace( np.matmul( Dr,Nop  ) )
    Nd = np.trace( np.matmul( Dd,Nop  ) )

    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )
    Sph = -np.trace( np.matmul(Dphf,aux) )

    El = np.trace( np.matmul( Dl,Hs ) )
    Er = np.trace( np.matmul( Dr,Hs ) )
    Ed = np.trace( np.matmul(Dd,Hs))

    return Nl.real, Ql.real, Qr.real, Qd.real, Qphlr.real, Sl.real, Sr.real,Sd.real, Sph.real, El.real, Er.real, Ed.real


E = 0
Ed = 0.
U0 = 0.
Uf = 50
g0 = 30/1000

eV = 6.5
mul1 = eV/2
mur1 = -eV/2
mud1 = 2

betaph = 1/100
betar,betad,betal = 1/100,1/2,1/100
gr,gd,gl = 1/100*(1/6),1/50,1/100
glU,grU = 1/100*(1/6),1/100
#gr,gd,gl = 1/100,0,1/100

#Ll = Dp(E,g0,U0,Uf,mul1,betal,gl) +  Dm(E,g0,U0,Uf,mul1,betal,gl)
#Lr = Dp(E,g0,U0,Uf,mur1,betar,gr) +  Dm(E,g0,U0,Uf,mur1,betar,gr)
#Ld = Dd(Ed,U0,mud1,betad,gd)

El =Er = E
J00 = 0.5
omegac0 = 1E-2
Ls = Dissipator(El,Er,g0,Ed,U0,Uf,mul1,mur1,mud1,betal,betar,betad,betaph,gl,glU,gr,grU,gd,J00,omegac0)
H = Hamiltonian(El,Er,Ed,U0,Uf,g0)

superop = Liouvillian(H,Ls)

#revisar la base
#|1,1,1>,|1,1,0>,|1,0,1>,|1,0,0>,|0,1,1>,|0,1,0>,|0,0,1>,|0,0,0>
alp,alp2,alp3,alp4,alp5 = 0.,0.,0.0,0.0,0.
a,b,c,d = 1j*alp,1j*alp2,1j*alp3,1j*alp4

rho0 = np.array([[1/8,0,0,0,0,0,0,0],
                 [0,1/8,a,0,d,0,0,0],
                 [0,-a,1/8,0,0,0,0,0],
                 [0,0,0,1/8,0,0,b,0],
                 [0,-d,0,0,1/8,0,0,0],
                 [0,0,0,0,0,1/8,c,0],
                 [0,0,0,-b,0,-c,1/8,0],
                 [0,0,0,0,0,0,0,1/8]])

rho1 = np.array([[0,0,0,0,0,0,0,0],
                 [0,1/4,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0],
                 [0,0,0,1/4,0,0,0,0],
                 [0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,1/4,0,0],
                 [0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,1/4]])






#for Uff in Ufs:
  #  calc = Propagate(rho0,superop,4000)   
 #   coheUf.append( abs(calc[5,3]))
                                           
#for Uu in Us:
  #  calc = Propagate(rho0,superop,4000)   
 #   coheU.append( abs(calc[5,3]))

#plt.plot(Ufs,coheUf)
#plt.show()


#plt.plot(Us,coheU)
#plt.show()

Ql = []
Qr = []
Qd = []

 
#Num = 1800
#J0s = np.logspace(-8,-3,Num)

Num = 3000
J0s = np.logspace(-8,1,Num)

Sls = []
Srs = []
Sds = []
Slr = []
Isl = []
Id = []
Els = []
Ers = []
Eds = []
Erl = []
cohev = []
concuv = []
Nls = []
Qphlist = []
Sphlist = []
Qlr = []
Jof = []

Probnt10 = []
Probnt20 = []
Probnt30 = []
Probnt40 = []
Probnt50 = []
Probnt60 = []
Probnt70 = []
Probnt80 = []
Imalphg = []
Imbetg = []

Ile = []
Ire = []
Iphs = []
Work = []
eff = []
effph = []
Realphg = []
Rebetg = []

#gf = 5/1000
g0fm = 1000
for J0 in J0s:
    mud0 = 2
    U00 = 40
    ev = 100
    #mud0 = 1-U00/2
    eps = 0.5
    Ed0f = mud0 - (1-eps)*U00    
    Uf0 = 500
    Elf = Erf = 0
    Delta,g,theta =energy(Elf,Erf,g0fm)
    Eup = (Elf+Erf)/2 + np.sqrt(Delta**2 + g**2)
    Emin = (Elf+Erf)/2 - np.sqrt(Delta**2 + g**2)
    dup,dupdag,dmin,dmindag = twolevel(Elf,Erf,g)
    omegac = 1E-2
    Ls0 = Dissipator(Elf,Erf,g,Ed0f,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,betaph,gl,glU,gr,grU,gd,J0,omegac)
    H0 = Hamiltonian(Elf,Erf,Ed0f,U00,Uf0,g)
    superop0 = Liouvillian(H0,Ls0)
    cal1f = Propagate(rho0,superop0,4000) 
    #Ll0 = Dp(E,g0,U00,Uf0,ev/2,betal,gl) + Dm(E,g0,U00,Uf0,ev/2,betal,gl)
    #Lr0 = Dp(E,g0,U00,Uf0,-ev/2,betar,gr) + Dm(E,g0,U00,Uf0,-ev/2,betar,gr)
    #Ld0 = Dd(Ed0,U00,mud0,betad,gd)
    #Hs0 = Hamiltonian(E,g0,Ed0,U00,Uf0)
    DLP,DRP = Dp(Eup,g,U00,Uf0,ev/2,betal,-ev/2,betar,gl,glU,gr,grU,dup,dupdag,dmin,dmindag,theta)
    DLM,DRM = Dm(Emin,g,U00,Uf0,ev/2,betal,-ev/2,betar,gl,glU,gr,grU,dmin,dmindag,dup,dupdag,theta)
    
    Laux = []
    Raux = []
    for l in DLP:
        Laux.append(l)
    for l in DLM:   
        Laux.append(l) 

    for r in DRP:
        Raux.append(r)
    for r in DRM:   
        Raux.append(r)

    DD = Dd(Elf,Erf,g,Ed0f,U00,mud0,betad,gd)
    Dph0 = Dphonon(Elf,Erf,g,betaph,J0,omegac)
    Ll0 = Laux
    Lr0 = Raux
    Ld0 = DD
    Nl0,Ql0,Qr0,Qd0,Qph0,Sl0,Sr0,Sd0,Sph0,El0,Er0,Ed0 = currents(Elf,Erf,g,H0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,Dph0,superop0,rho0,30000)
    
    #hay que medir de otra forma la coherencia 
    #elegir bien como medir coherencia y entrelazamient
    #alp0 = np.trace( np.matmul( np.matmul(dldag,dr),cal1f) )
    alp0 = v100.conjugate().T @ cal1f @ v010
    alp1 = v101.conjugate().T @ cal1f @ v011
    #print(alp0[0,0],alp1[0,0])
    P0 = v00.conjugate().T @ cal1f @ v00 + v001.conjugate().T @ cal1f @ v001
    PD = v110.conjugate().T @ cal1f @ v110 + v111.conjugate().T @ cal1f @ v111
    cohev.append(abs(alp0[0,0]) + abs(alp1[0,0]) )
    cohesum = abs(alp0[0,0] + alp1[0,0])
    concurrencef = 2*cohesum - 2*np.sqrt(abs(P0[0,0]*PD[0,0]))
    if (concurrencef > 0):
        concuv.append(concurrencef)
    else:
        concuv.append(0)
    Nls.append(Nl0)
    Ql.append(Ql0)
    Qr.append(Qr0)
    Qd.append(Qd0)
    Qphlist.append(Qph0)
    Jof.append(J0/(betaph*gl))

    Qlr.append(Ql0 + Qr0 + Qph0)
    sigmal = Sl0 - betal*Ql0
    sigmar = Sr0 - betar*Qr0
    sigmaph = Sph0 - betaph*Qph0
    
    Probnt10.append(cal1f[0,0].real )
    Probnt20.append(cal1f[1,1].real )
    Probnt30.append(cal1f[2,2].real )
    Probnt40.append(cal1f[3,3].real )
    Probnt50.append(cal1f[4,4].real )
    Probnt60.append(cal1f[5,5].real )
    Probnt70.append(cal1f[6,6].real )
    Probnt80.append(cal1f[7,7].real )
    Imalphg.append(2*g0fm*alp0[0,0].real)
    Imbetg.append(2*g0fm*alp1[0,0].real)

    Realphg.append(alp0[0,0].real)
    Rebetg.append(alp1[0,0].real)

    Sls.append((Sl0 - betal*Ql0))
    Srs.append((Sr0 - betar*Qr0))
    Sds.append(Sd0 - betad*Qd0)
    Slr.append( sigmal + sigmar )
    Sphlist.append( sigmaph )
    Isl.append(-Sl0 - Sr0)
    Id.append(-Sd0)
    Els.append(El0)
    Ers.append(Er0)
    Eds.append(Ed0)
    Erl.append(El0 + Er0)
    Ile.append(-Sl0)
    Ire.append(-Sr0)
    Iphs.append(-Sph0)
    Wlr0 = ev*Nl0
    Work.append(Wlr0)
    eff.append( abs(Wlr0)/(Ql0+Qr0) )
    effph.append( abs(Wlr0)/(Ql0+Qr0+Qph0) )

    print(J0)


fig, ax10 = plt.subplots(
    sharex=True,
    figsize=(3.39, 3.2)
)

LINE_W = 1.6
LABEL_FS = 9
TICK_FS = 8
PANEL_FS = 9

# ---------- Panel (a)
ax10.plot(Jof, Realphg, color='blue',   lw=LINE_W, label=r'$2gRe(a)$')
ax10.plot(Jof, Rebetg, color='orange', lw=LINE_W, label=r'$2gRe(b)$')

#ax10.plot(Jof5, Imalphg5, color='brown',  lw=LINE_W, label=r'$1$')
ax10.set_xlabel(r'$g/\kappa_L$', fontsize=LABEL_FS)
ax10.set_ylabel(r'$\mathcal{C}_{l_{1}}$', fontsize=LABEL_FS)
ax10.set_xscale('log')
ax10.tick_params(direction='in', which='both', labelsize=TICK_FS)
#ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes,
          #fontsize=PANEL_FS, fontweight='bold')

ax10.legend(
    fontsize=6.5,
    frameon=True,
    ncol=1,
    loc='upper right',
    bbox_to_anchor=(1, 0.99)
)



# ---------- Spines

for spine in ax10.spines.values():
    spine.set_linewidth(0.8)

plt.tight_layout(pad=0.4)
plt.show()


###ideas: calcular <100|\rho|010> y <101|\rho|011> 
#np.savez("phonong=0.1b100.npz", Jof=Jof, Id=Id,Ile =Ile,Ire = Ire, Iphs = Iphs,cohes=cohev, concv = concuv, Nls = Nls,Work=Work, eff=eff,effph=effph)
np.savez("phonong=1000b100sec.npz", Jof=Jof, Id=Id,Ile =Ile,Ire = Ire, Iphs = Iphs,cohes=cohev, concv = concuv, Nls = Nls,Work=Work, eff=eff,effph=effph)

np.savez("phonong=1000probb100sec.npz", Jof=Jof, Probnt10=Probnt10,Probnt20=Probnt20,Probnt30=Probnt30,Probnt40=Probnt40,Probnt50=Probnt50,Probnt60=Probnt60,Probnt70=Probnt70,Probnt80=Probnt80, Imalphg=Imalphg, Imbetg=Imbetg)