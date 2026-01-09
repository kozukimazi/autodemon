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

def Dissipator(El,Er,g,Ed,U,Uf,mul,mur,mud,betal,betar,betad,gl,glU,gr,grU,gammad):
    Delta, g, theta = energy(El,Er,g)
    eps = (El+Er)/2
    Eup, Emin = eps + np.sqrt(Delta**2 + g**2), eps - np.sqrt(Delta**2 + g**2)
    
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    DLP,DRP = Dp(Eup,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dup,dupdag,dmin,dmindag,theta)
    DLM,DRM = Dm(Emin,g,U,Uf,mul,betal,mur,betar,gl,glU,gr,grU,dmin,dmindag,dup,dupdag,theta)
   
    DD = Dd(El,Er,g,Ed,U,mud,betad,gammad)

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
    a1 = Eup*np0 + Emin*nm + Ed*nd
    #a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(np0,nd) +  np.matmul(nm,nd) ) + Uf*np.matmul(np0,nm) 
    return a1+a3

def currents(El,Er,g,Hs,mul,mur,mud,Ll,Lr,Ld,superop,rho0,t):
    dup,dupdag,dmin,dmindag = twolevel(El,Er,g)
    nm = np.matmul(dmindag,dmin)
    np0 = np.matmul(dupdag,dup)
    Nop = np0 + nm + nd
    rhof = Propagate(rho0,superop,t)

    Dl = Dissipate(Hs,Ll,rhof)
    Dr = Dissipate(Hs,Lr,rhof)
    Dd = Dissipate(Hs,Ld,rhof)

    Qopl = (Hs - mul*Nop)
    Qopr = (Hs - mur*Nop)
    Qopd = (Hs - mud*Nop)
    aux = logm(rhof)
    Ql = np.trace( np.matmul( Dl,Qopl  ) )
    Qr = np.trace( np.matmul( Dr,Qopr  ) )
    Qd = np.trace( np.matmul( Dd,Qopd  ) )

    Nl = np.trace( np.matmul( Dl,Nop  ) )
    Nr = np.trace( np.matmul( Dr,Nop  ) )
    Nd = np.trace( np.matmul( Dd,Nop  ) )

    Wl = mul*np.trace( np.matmul( Dl,Nop  ) )
    Wr = mur*np.trace( np.matmul( Dr,Nop  ) )
    


    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )

    El = np.trace( np.matmul( Dl,Hs  ) )
    Er = np.trace( np.matmul( Dr,Hs  ) )
    Ed = np.trace( np.matmul(Dd,Hs))

    return Nl.real, Wl.real,Wr.real,Ql.real, Qr.real, Qd.real, Sl.real, Sr.real,Sd.real, El.real, Er.real, Ed.real

betar,betad,betal = 1/100,1/2,1/100
gr,gd,gl = 1/100*(1/6),1/50,1/100
glU,grU = 1/100*(1/6),1/100

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


cohev = []
concuv = []
Ilr = []
Jlr = []
Wlr = []
Nlr = []
g0m = []

#Num = 80000
#Num = 10
#g0s = np.linspace(0,1,Num)
Num = 5000
g0s = np.logspace(-6,3,Num)
for g0f in g0s:
    #aqui hay funciones que no deberian estar para optimizarlo
    ev = 100
    print(g0f)
    mud0 = 2
    U00 = 40
    #mud0 = 1-U00/2
    eps = 0.5
    Ed0f = mud0 - (1-eps)*U00    
    Uf0 = 500
    Elf = Erf = 0
    Delta,g,theta = energy(Elf,Erf,g0f)
    Eup = (Elf+Erf)/2 + np.sqrt(Delta**2 + g0f**2)
    Emin = (Elf+Erf)/2 - np.sqrt(Delta**2 + g0f**2)
    dup,dupdag,dmin,dmindag = twolevel(Elf,Erf,g0f)
    
    Ls0 = Dissipator(Elf,Erf,g0f,Ed0f,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,gl,glU,gr,grU,gd)
    H0 = Hamiltonian(Elf,Erf,Ed0f,U00,Uf0,g0f)
    superop0 = Liouvillian(H0,Ls0)
    cal1f = Propagate(rho0,superop0,10000) 
    #Ll0 = Dp(E,g0,U00,Uf0,ev/2,betal,gl) + Dm(E,g0,U00,Uf0,ev/2,betal,gl)
    #Lr0 = Dp(E,g0,U00,Uf0,-ev/2,betar,gr) + Dm(E,g0,U00,Uf0,-ev/2,betar,gr)
    #Ld0 = Dd(Ed0,U00,mud0,betad,gd)
    #Hs0 = Hamiltonian(E,g0,Ed0,U00,Uf0)
    DLP,DRP = Dp(Eup,g0f,U00,Uf0,ev/2,betal,-ev/2,betar,gl,glU,gr,grU,dup,dupdag,dmin,dmindag,theta)
    DLM,DRM = Dm(Emin,g0f,U00,Uf0,ev/2,betal,-ev/2,betar,gl,glU,gr,grU,dmin,dmindag,dup,dupdag,theta)
    
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

    DD = Dd(Elf,Erf,g0f,Ed0f,U00,mud0,betad,gd)
    Ll0 = Laux
    Lr0 = Raux
    Ld0 = DD
    
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
    g0m.append(g0f/gl)    

    Nl0,Wl0,Wr0,Ql0,Qr0,Qd0,Sl0,Sr0,Sd0,El0,Er0,Ed0 = currents(Elf,Erf,g,H0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,superop0,rho0,30000)
    Ilr.append(-Sl0-Sr0)
    Jlr.append(Ql0+Qr0)
    Wlr.append(Wl0+Wr0)

archivo = open("globalcouplarge","w")
decimal_places = 7
total_width = 8
format_str = f"{{:.{decimal_places}f}}" 
#format_str = f"{{:{total_width}.{decimal_places}f}}"
for i in range(Num):
    archivo.write( format_str.format(g0m[i])) #guarda el grado del nodo
    #archivo.write(str(xs[i])) 
    archivo.write(" ") 
    archivo.write( format_str.format(concuv[i]))
    archivo.write(" ") 
    archivo.write( format_str.format(cohev[i]))
    archivo.write(" ")
    archivo.write( format_str.format(Jlr[i]))
    archivo.write(" ")
    archivo.write( format_str.format(Wlr[i]))
    archivo.write(" ")
    archivo.write( format_str.format(Ilr[i]))
    archivo.write("\n")
