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

def Logaritm(p10i,p01i,cohei):
     
    lamb1p = (p10i+p01i)/2 + (np.sqrt( (p10i-p01i)**2 + 4*(abs(cohei)**2) ))/2
    lamb1m = (p10i+p01i)/2 - (np.sqrt( (p10i-p01i)**2 + 4*(abs(cohei)**2) ))/2
    delt0 = (p10i-p01i)/2
    costh = abs(cohei)/np.sqrt( delt0**2 + abs(cohei)**2 )
    ##cos(\theta/2)
    costh2 = np.sqrt( (1+ costh)/2 )
    fase = cohei/abs(cohei)
    a = (costh2**2)*np.log(lamb1p) + (1-(costh2**2))*np.log(lamb1m) 
    b = (1-(costh2**2))*np.log(lamb1p) + (costh2**2)*np.log(lamb1m)
    ci = (np.sqrt(1-costh**2)/2)*fase*(np.log(lamb1p) - np.log(lamb1m) )

    return a,b,ci
    
def currents(Htd,mul,mur,mud,Ll,Lr,Ld,superop,rho0,t):
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

    Qopl = (Htd - mul*Nop)
    Qopr = (Htd - mur*Nop)
    Qopd = (Htd - mud*Nop)

    aux = logm(rhof)
    Ql = np.trace( np.matmul( Dl,Qopl  ) )
    Qr = np.trace( np.matmul( Dr,Qopr  ) )
    Qd = np.trace( np.matmul( Dd,Qopd  ) )

    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )

    El = np.trace( np.matmul( Dl,Htd  ) )
    Er = np.trace( np.matmul( Dr,Htd  ) )
    Ed = np.trace( np.matmul(Dd,Htd))

    Wl = mul*np.trace(np.matmul( Dl, Nop ))
    Wr = mur*np.trace(np.matmul( Dr, Nop ))
    Wd = mud*np.trace(np.matmul( Dd, Nop ))

    Nl = np.trace( np.matmul(Dl,Nop ) )
    Nr = np.trace(np.matmul(Dr,Nop ))
    Nd = np.trace(np.matmul( Dd,Nop ))
    return Ql.real, Qr.real, Qd.real, Sl.real, Sr.real,Sd.real, El.real, Er.real, Ed.real, Wl.real,Wr.real,Wd.real,cohe,concf,Nl.real,Nr.real,Nd.real


E = 0
U0 = 40.
Uf = 500
g0 = 5/1000
#g0 =1/1000 puede ser mejor

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

gr,grU = (1/100)*(1/6), 1/100
gl,glU = 1/100, (1/100)*(1/6)
gd,gdU = 1/50,1/50 
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

Ll = Dl(E,U0,Uf,mul1,betal,gl,glU)
Lr = Dr(E,U0,Uf,mur1,betar,gr,grU)
Ld = Dd(Ed,U0,mud1,betad,gd,gdU)



Ls = Dissipator(E,Ed,U0,Uf,mul1,mur1,mud1,betal,betar,betad,gl,glU,gr,grU,gd,gdU)
H = Hamiltonian(E,Ed,U0,Uf,g0)

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



times = np.linspace(0,2600,2000)
Probnt1 = []
Probnt2 = []
Probnt3 = []
Probnt4 = []
Probnt5 = []
Probnt6 = []
Probnt7 = []
Probnt8 = []
traza = []
cohe = []
concu = []
for ts in times:
    cal1 = Propagate(rho0,superop,ts)
    auxp = np.matmul(dldag,dr)
    alp = np.trace( np.matmul(auxp,cal1) )
    tot = np.trace(cal1)
    
    traza.append(tot)
    Probnt1.append(cal1[0,0].real )
    Probnt2.append(cal1[1,1].real )
    Probnt3.append(cal1[2,2].real )
    Probnt4.append(cal1[3,3].real )
    Probnt5.append(cal1[4,4].real )
    Probnt6.append(cal1[5,5].real )
    Probnt7.append(cal1[6,6].real ) 
    Probnt8.append(cal1[7,7].real ) 
    cohe.append(abs(cal1[5,3]) + abs(cal1[4,2]) )
    cohesum = abs(cal1[5,3] + cal1[4,2])
    PD = cal1[0,0].real + cal1[1,1].real 
    P0 = cal1[7,7].real + cal1[6,6].real 
    concurrence = 2*cohesum - 2*np.sqrt(P0*PD) 
    concu.append(concurrence)

    
    #print(cal1[4,2])
    #print(cal1[2,4])
    #cal2 = DistR(H,Ls,Lr,Ll,rho0,ts,-2)
    #Probnt2.append(cal2.real)

plt.plot(times,Probnt1,label = r'$\rho_{111}$')
plt.plot(times,Probnt2, label = r'$\rho_{110}$')
plt.plot(times,Probnt3, label = r'$\rho_{101}$')
plt.plot(times,Probnt4, label = r'$\rho_{100}$')
plt.plot(times,Probnt5, label = r'$\rho_{011}$')
plt.plot(times,Probnt6, label = r'$\rho_{010}$')
plt.plot(times,Probnt7, label = r'$\rho_{001}$')
plt.plot(times,Probnt8, label = r'$\rho_{000}$')
plt.legend(loc = "upper right",fontsize=15)
plt.xlabel(r'$t$',fontsize=25)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.show()


plt.plot(times,cohe)
plt.xlabel(r'$t$', fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$', fontsize = 20)
plt.show()
#plt.scatter(times,Probnt2,label = "Baño_d")
plt.plot(times,concu)
plt.xlabel(r'$t$', fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$', fontsize = 20)
plt.show()

values =  Propagate(rho0,superop,3000)

plt.imshow(values.imag)
plt.colorbar()
plt.show()
    
Ufs = np.linspace(7,40,100)
Us = np.linspace(1,6,50)
Num = 80
eVs0 = np.linspace(0,1800,Num)


Probnt1f = []
Probnt2f = []
Probnt3f = []
Probnt4f = []
Probnt5f = []
Probnt6f = []
Probnt7f = []
Probnt8f = []
trazaf = []
cohef = []
concuf = []
eVs = []
a0s = []
b0s = []
c0s = []
a1s = []
b1s = []
c1s = []
Realp = []
Rebet = []
infoqm = []
imalp = []
imbet = []
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
    rhof = Propagate(rho0,superop0,40000)
    trazaf.append(tot)

    p100,p010,alp = rhof[3,3].real,rhof[5,5].real,rhof[3,5]
    p101,p011,bet = rhof[2,2].real,rhof[4,4].real,rhof[2,4]
    a0,b0,c0 = Logaritm(p100,p010,alp)
    a1,b1,c1 = Logaritm(p101,p011,bet)
    Realp.append((alp*c0.conj()).real )
    Rebet.append((bet.conj()*c1).real )
    infoqm.append(grU*(1-fermi(0,-ev/2,betar))*(bet.conj()*c1).real )
    Probnt1f.append(rhof[0,0].real )
    Probnt2f.append(rhof[1,1].real )
    Probnt3f.append(rhof[2,2].real )
    Probnt4f.append(rhof[3,3].real )
    Probnt5f.append(rhof[4,4].real )
    Probnt6f.append(rhof[5,5].real )
    Probnt7f.append(rhof[6,6].real ) 
    Probnt8f.append(rhof[7,7].real ) 
    imalp.append(rhof[3,5].imag)
    imbet.append(rhof[2,4].imag)
    cohef.append(abs(rhof[5,3]) + abs(rhof[4,2]) )
    cohesumf = abs(rhof[5,3] + rhof[4,2])
    PDf = rhof[0,0].real + rhof[1,1].real 
    P0f = rhof[7,7].real + rhof[6,6].real 
    concurrencef = 2*cohesumf - 2*np.sqrt(P0f*PDf) 
    concuf.append(concurrencef)
    eVs.append(ev*betal)
    a0s.append(a0)
    b0s.append(b0)
    c0s.append(abs(c0))
    a1s.append(a1)
    b1s.append(b1)
    c1s.append(abs(c1))

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


plt.plot(eVs,cohef)
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{l_{1}}$', fontsize = 20)
plt.show()
#plt.scatter(times,Probnt2,label = "Baño_d")
plt.plot(eVs,concuf)
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.ylabel(r'$\mathcal{C}_{on}$', fontsize = 20)
plt.show()



plt.plot(eVs,a0s,label=r'$a_{0}$')
plt.plot(eVs,b0s,label = r'$b_{0}$')
plt.plot(eVs,c0s,label = r'$c_{0}$')
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.legend()
plt.show()
#plt.scatter(times,Probnt2,label = "Baño_d")
plt.plot(eVs,a1s,label=r'$a_{1}$')
plt.plot(eVs,b1s,label = r'$b_{1}$')
plt.plot(eVs,c1s,label = r'$c_{1}$')
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.legend()
plt.show()


plt.plot(eVs,Realp,label=r'$Re(\alpha c^{*}_{0})$')
plt.plot(eVs,Rebet,label = r'$Re(\beta^{*}c_{1})$')
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.legend()
plt.show()

plt.plot(eVs,infoqm,label=r'$I_{qm}$')
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.legend()
plt.show()

plt.plot(eVs,imalp,lw = 4,label=r'$\text{Im}(\alpha)$')
plt.plot(eVs,imbet,lw = 4,label=r'$\text{Im}(\beta)$')
plt.xlabel(r'$eV/T$', fontsize = 20)
plt.xticks(fontsize=17)  
plt.yticks(fontsize=17)
plt.legend(fontsize = 15)
plt.show()