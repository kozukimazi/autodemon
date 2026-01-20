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

tot = anticonmutador(dddag,dd)
totl = anticonmutador(dldag,dl)
totr = anticonmutador(drdag,dr)


nd = np.matmul(dddag,dd)
nl = np.matmul(dldag,dl)
nr = np.matmul(drdag,dr)


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

def Dphonon(betaph,J0,omegac):
    gam = gammaph(0,omegac,J0,betaph)
    aux = np.sqrt(gam) *(np.matmul(dldag, dr) + np.matmul(drdag, dl))
    return [aux]


def Dissipator(E,Ed,U,Uf,mul,mur,mud,betal,betar,betad,betaph,gammal,gammalU,gammar,gammarU,gammad,gammadU,J0,omegac):
    DR = Dr(E,U,Uf,mur,betar,gammar,gammarU)
    DL = Dl(E,U,Uf,mul,betal,gammal,gammalU)
    DD = Dd(Ed,U,mud,betad,gammad,gammadU)
    Dph = Dphonon(betaph,J0,omegac)
    tot = []
    for l in DL:
        tot.append(l)
    for r in DR:
        tot.append(r)
    for d in DD:
        tot.append(d)    

    for p in Dph:
        tot.append(p)
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

def adjLiouvillian( H,Ls,O,hbar = 1):
    
    cohe = (1j/hbar)* (np.matmul(H,O ) - np.matmul(O,H))
    dissip = sum([ np.matmul( L.conjugate().T,np.matmul(O,L) ) - 1/2*( anticonmutador( np.matmul(L.conjugate().T,L),O ) )  for L in Ls] )
    return cohe + dissip

def adjLiouvillianT( H,Ls,hbar = 1):
    d = len(H)
    cohe = (1j/hbar)* (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d)))
    dissip = sum([ np.kron( L.transpose(),L.transpose().conjugate(), ) - 1/2*( np.kron((L.conjugate().transpose().dot(L)).T,np.eye(d)) + np.kron(np.eye(d),(L.conjugate().transpose().dot(L))) )  for L in Ls] )
    return cohe + dissip    


def Propagate(rho0,superop,t):
    d = len(rho0)
    propagator = expm (superop *t)
    vec_rho_t = propagator @ np.reshape(rho0,(d**2,1))
    return np.reshape( vec_rho_t, (d,d) )

def Hamiltonian(El,Er,Ed,U,Uf,g):
    a1 = El*nl + Er*nr + Ed*nd
    a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(nl,nd) +  np.matmul(nr,nd) ) + Uf*np.matmul(nl,nr) 
    return a1+a2+a3

def Htd(El,Er,Ed,U,Uf):
    Htdl = El*nl + Er*nr + Ed*nd + U*(np.matmul(nr,nd) + np.matmul(nl,nd)) + Uf*np.matmul(nr,nl) 
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


def currents(Htd,g0,mul,mur,mud,Ll,Lr,Ld,Dph,superop,rho0,t):
    Nop = nl + nr + nd
    Iquantum = -1j*g0*(np.matmul(dldag,dr) - np.matmul(drdag,dl))
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
    Dph0 = Dissipate(Htd,Dph,rhof)

    Qopl = (Htd - mul*Nop)
    Qopr = (Htd - mur*Nop)
    Qopd = (Htd - mud*Nop)
    A = np.matmul(dldag,dr) + np.matmul(drdag,dl)

    aux = logm(rhof)
    Ql = np.trace( np.matmul( Dl,Qopl  ) )
    Qr = np.trace( np.matmul( Dr,Qopr  ) )
    Qd = np.trace( np.matmul( Dd,Qopd  ) )
    Qph = np.trace( np.matmul( Dph0,Htd  ) )

    Nqm = np.trace( np.matmul( Iquantum, rhof ) )

    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )
    Sph = -np.trace( np.matmul(Dph0,aux) )

    El = np.trace( np.matmul( Dl,Htd  ) )
    Er = np.trace( np.matmul( Dr,Htd  ) )
    Ed = np.trace( np.matmul(Dd,Htd))
    
    Wl = mul*np.trace(np.matmul( Dl, Nop ))
    Wr = mur*np.trace(np.matmul( Dr, Nop ))
    Wd = mud*np.trace(np.matmul( Dd, Nop ))
    Wlr = Wl + Wr
    Nl = np.trace( np.matmul(Dl,Nop ) )
    Nr = np.trace(np.matmul(Dr,Nop ))
    Nd = np.trace(np.matmul( Dd,Nop ))
    Act = np.trace(np.matmul(A,np.matmul(rhof,A)))
    return Ql.real, Qr.real, Qd.real, Qph.real, Sl.real, Sr.real,Sd.real, Sph.real, El.real, Er.real, Ed.real, Wl.real,Wr.real,Wd.real,cohe,concf,Nl.real,Nr.real,Nd.real,Act.real,Nqm.real,Wlr.real

###################################
#############parametros##############
###################################
E = 0
U0 = 40.
Uf = 500
#caso solo fonones
g0 =1/1000
#g0 = 5/1000
#g0 = 600
#g0 = 1/1000, pasa algo muy interesante con el entrelazamiento, muere y reaparece
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
betaph = 1/100
J0, omegac = 0.000005, 1E-2 
gr,grU = (1/100)*(1/6), 1/100
gl,glU = 1/100, (1/100)*(1/6)
gd,gdU = 1/50,1/50
#gd,gdU = 1/400,1/400  
#gr,grU = (1/100), 1/100
#gl,glU = 1/100, (1/100)
#gd,gdU = 1/100,1/100 
#al crecer gd,gdu mayor a las gls y grs
# se aumenta el rango de eV en el que hay rompimiento aparente 
#segunda ley
#al disminuir gd,gdu se disminuye el rango de eV

#
#gr,gd,gl = 1/100,1/300,1/100
#gr,gd,gl = 1/100,0,1/100

El=Er = 0


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

Num = 3000
J0s = np.logspace(-8,1,Num)

#Num = 1800
#J0s = np.logspace(-8,-3,Num)

Isl = []
Id = []
Iphs = []
Ile = []
Ire = []
cohes = []
concv = []
Ilrnew = []
Jof = []
Nls = []
Acts = []
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
Nqms = []
Nltotal = []
Nlqm = []
Work = []
Qls = []
Qrs = []
Qds = []
eff = []
effph = []
for J0 in J0s:
    mud0 = 2
    U00 = 40 #10
    #mud0 = 1-U00/2
    #Con Ed0 = mud0 -U00/2,E0=4 hay flujo de energia pero un orden menor al de
    #flujo de informacion
    #Ed0 = 1
    eps = 0.5
    Ed0 = mud0 - (1-eps)*U00
    ev = 100
    #Ed0 = mud0 -U00/2
    Uf0 = 500 #50
    #Probar condicion (U00/E0)<<1,Strasberg
    El0 = Er0 = E0 = 0
    Ls0 = Dissipator(E0,Ed0,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,betaph,gl,glU,gr,grU,gd,gdU,J0,omegac)
    Ls0p = Dissipator(E0,Ed0,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,betaph,0,0,0,0,0,0,J0,omegac)
    H0 = Hamiltonian(El0,Er0,Ed0,U00,Uf0,g0)
    H0f = Hamiltonian(El0,Er0,Ed0,U00,Uf0,0)
    superop0 = Liouvillian(H0,Ls0)
    superop0adj = adjLiouvillian(H0,Ls0p,nl)
    superop1adj = adjLiouvillian(H0,0*Ls0p,nl)
    superop0adj2 = adjLiouvillian(H0,Ls0,nr)
    final = superop0adj + superop0adj2
    adLin = adjLiouvillianT(H0,Ls0)
    #nlt = Propagate(nl,adLin,40000)
    Ll0 = Dl(E0,U00,Uf0,ev/2,betal,gl,glU)
    Lr0 = Dr(E0,U00,Uf0,-ev/2,betar,gr,grU)
    Ld0 = Dd(Ed0,U00,mud0,betad,gd,gdU)
    Lph0 = Dphonon(betaph,J0,omegac)
    Htd0 =  Htd(El0,Er0,Ed0,U00,Uf0)
    rhof = Propagate(rho0,superop0,40000)
    Nlt = np.trace( np.matmul(superop0adj,rhof ) )
    Nltqm = np.trace( np.matmul(superop1adj,rhof ) )
    #Nlt = np.trace( np.matmul(nlt,rho0 ) )
    Ql0,Qr0,Qd0,Qph0,Sl0,Sr0,Sd0,Sphf0,El0,Er0,Ed0,Wl0,Wr0,Wd0,cohe0,concu0,Nl0,Nr0,Nd0,Act0,Nq0,Wlr0 = currents(Htd0,g0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,Lph0,superop0,rho0,40000)
    #cohev.append(abs(rhof[5,3]) + abs(rhof[4,2]) )
    concv.append(concu0)    
    sigmal = Sl0 - betal*Ql0
    sigmar = Sr0 - betar*Qr0

    Isl0 = -Sl0 - Sr0
    Id0 = -Sd0
    Il0 = -Sl0
    Ir0 = -Sr0
    Iph0 = -Sphf0
    Isl.append(Isl0)
    Id.append(-Sd0)
    Iphs.append(Iph0)
    Ilrnew.append(Il0 + Ir0 + Iph0)
    Ile.append(Il0)
    Ire.append(Ir0)
    cohes.append(cohe0)
    Nls.append(Nl0)
    Jof.append(J0/(betaph*gl))
    Acts.append(J0*Act0)
    Probnt10.append(rhof[0,0].real )
    Probnt20.append(rhof[1,1].real )
    Probnt30.append(rhof[2,2].real )
    Probnt40.append(rhof[3,3].real )
    Probnt50.append(rhof[4,4].real )
    Probnt60.append(rhof[5,5].real )
    Probnt70.append(rhof[6,6].real )
    Probnt80.append(rhof[7,7].real )
    Imalphg.append(2*g0*rhof[5,3].imag)
    Imbetg.append(2*g0*rhof[4,2].imag)
    Nqms.append(Nq0)
    Nltotal.append(Nlt)
    Nlqm.append(Nltqm)
    Work.append(Wlr0)
    Qls.append(Ql0)
    Qrs.append(Qr0) 
    Qds.append(Qd0)
    eff.append( abs(Wlr0)/(Ql0+Qr0) )
    effph.append( abs(Wlr0)/(Ql0+Qr0+Qph0) )
    print(J0)


plt.plot(Jof,Id, color='red',lw=3, label = r'$\dot{I}_{D}$')
plt.plot(Jof,Iphs, color='orange',lw=3, label = r'$\dot{I}_{Ph}$')
plt.plot(Jof,Ile, color='black',lw=3, label = r'$\dot{I}_{Le}$')
plt.plot(Jof,Ire, color='blue',lw=3, label = r'$\dot{I}_{Re}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.ylabel(r'$\dot{I}_{i}$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper left")
plt.xscale("log")
plt.show()

plt.plot(Jof,concv, color='green',lw=3, label = 'Concurrence')
plt.plot(Jof,cohes, color='purple',lw=3, label = 'Coherence')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof,Nls, color='green',lw=3, label = r'$\dot{N}_{L}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof,Nlqm, color='green',lw=3, label = r'$\hat{I}$')
plt.plot(Jof,Nltotal, color='red',lw=3,linestyle = '--', label = r'$\hat{I} + \hat{I}_{diss}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


plt.plot(Jof,Qls, color='green',lw=3, label = r'$\dot{Q}_{L}$')
plt.plot(Jof,Qrs, color='red',lw=3,linestyle = '--', label = r'$\dot{Q}_{R}$')
plt.plot(Jof,Qds, color='blue',lw=3, label = r'$\dot{Q}_{D}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()

plt.plot(Jof,eff, color='green',lw=3, label = r'$\eta$')
plt.plot(Jof,effph, color='red',lw=3, label = r'$\eta_{ph}$')
plt.xlabel(r'$J_{0}/(\beta_{ph}\gamma_{L})$',fontsize = 20)
plt.xticks(fontsize=17)  # X-axis tick labels
plt.yticks(fontsize=17)
plt.legend(fontsize=15,loc = "upper right")
plt.xscale("log")
plt.show()


# Create subplots (1 row, 2 columns)
fig, (ax10, ax20) = plt.subplots(2, 1,sharex=True, figsize=(4, 9),constrained_layout=True)  # 1 row, 2 columns


#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
ax10.plot(Jof,Iphs, color='red',lw=3, label = r'$\dot{I}_{ph}$')
ax10.plot(Jof,Ile, color='black',lw=3, label = r'$\dot{I}_{Le}$')
ax10.plot(Jof,Ire, color='blue',lw=3, label = r'$\dot{I}_{Re}$')
ax10.set_ylabel(r'$\dot{I}_{i}$',fontsize = 22)
ax10.legend(fontsize=17,loc = "center left")
ax10.set_xscale('log')  
ax10.tick_params(labelbottom=False,labelsize = 18)
ax10.text(0.9, 0.93, '(a)', transform=ax10.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


ax20.plot(Jof,cohes,label = r'$\mathcal{C}_{l_{1}}$', color = 'b',lw = 3)
ax20.plot(Jof,concv, label = r'$\mathcal{C}_{on}$', color = 'r',lw=3)  
ax20.set_xlabel(r'$J_{0}/(\beta_{ph}\kappa_{L})$',fontsize = 20)   
ax20.set_xscale("log")
ax20.legend(fontsize=17, loc = "upper left") 
ax20.set_ylabel("Coherencia y entrelazamiento",fontsize = 22)
ax20.tick_params(labelsize=18)  # font size of tick labels 
ax20.text(0.9, 0.93, '(b)', transform=ax20.transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()  # Avoids overlapping labels
plt.show()



np.savez("phonong=10^{-3}b100ne.npz", Jof=Jof, Id=Id,Ile =Ile,Ire = Ire, Iphs = Iphs,cohes=cohes, concv = concv, Nls = Nls,Acts=Acts,Nlqm=Nlqm,Nltotal=Nltotal,Work=Work, eff=eff,effph=effph)

#np.savez("phonong=3_10^{-3}probb100.npz", Jof=Jof, Probnt10=Probnt10,Probnt20=Probnt20,Probnt30=Probnt30,Probnt40=Probnt40,Probnt50=Probnt50,Probnt60=Probnt60,Probnt70=Probnt70,Probnt80=Probnt80, Imalphg=Imalphg, Imbetg=Imbetg)
#volver a sacar 3_10^{-3}(16/12/25)
