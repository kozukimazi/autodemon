import numpy as np
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

def Logaritm(rho):
    po,p1,p2,p3 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p4,p5,p6,p7 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    bet,alp = rho[2,4],rho[3,5]   

    p00,p11,p66,p77 = np.log(po),np.log(p1),np.log(p6),np.log(p7)
    lamb1p = (p2+p3)/2 + (np.sqrt( (p2-p3)**2 + 4*(abs(alp)**2) ))/2
    lamb1m = (p2+p3)/2 - (np.sqrt( (p2-p3)**2 + 4*(abs(alp)**2) ))/2
    delt1 = (p2-p3)/2 
    lamb0p = (p4+p5)/2 + (np.sqrt( (p4-p5)**2 + 4*(abs(bet)**2) ))/2
    lamb0m = (p4+p5)/2 - (np.sqrt( (p4-p5)**2 + 4*(abs(bet)**2) ))/2
    delt0 = (p4-p5)/2
    #if (abs(delt)< 1E-10):
     #   th = np.pi/2
    #else:    
    th1=arctanaux(abs(alp)/delt1)
    fase1 = cmath.phase(alp)
    fase = np.exp(1j*fase1)
    th0=arctanaux(abs(bet)/delt0)
    fase0 = cmath.phase(bet)
    fase00 = np.exp(1j*fase0)
    #if (abs(alp) > 0 ):
     #   fase += alp/abs(alp)
    alp1 = (np.cos(th1/2)**2)*np.log(lamb1p) + (np.sin(th1/2)**2)*np.log(lamb1m)
    bet1 = (np.sin(th1/2)**2)*np.log(lamb1p) + (np.cos(th1/2)**2)*np.log(lamb1m)   
    c1 =  np.cos(th1/2)*np.sin(th1/2)*(fase)*(np.log(lamb1p) - np.log(lamb1m) )
    
    alp0 = (np.cos(th0/2)**2)*np.log(lamb0p) + (np.sin(th0/2)**2)*np.log(lamb0m)
    bet0 = (np.sin(th0/2)**2)*np.log(lamb0p) + (np.cos(th0/2)**2)*np.log(lamb0m)   
    c0 =  np.cos(th0/2)*np.sin(th0/2)*(fase00)*(np.log(lamb0p) - np.log(lamb0m) )

    matrix = np.zeros((8,8), dtype = np.complex_)
    matrix[0,0],matrix[1,1],matrix[6,6],matrix[7,7] = p00,p11,p66,p77
    matrix[2,2],matrix[3,3],matrix[4,4],matrix[5,5] = alp1,bet1,alp0,bet0
    matrix[2,4],matrix[4,2] = c1,c1.conj()
    matrix[3,5],matrix[5,3] = c0,c0.conj()

    return matrix 

def Log2(rho):
    p111,p110,p101,p100 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p011,p010,p001,p000 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    alp,bet = rho[2,4],rho[3,5]   

    lamb1p = (p101+p011)/2 + (np.sqrt( (p101-p011)**2 + 4*(abs(alp)**2) ))/2
    lamb1m = (p101+p011)/2 - (np.sqrt( (p101-p011)**2 + 4*(abs(alp)**2) ))/2
    delt1 = (p101-p011)/2 
    lamb0p = (p100+p010)/2 + (np.sqrt( (p100-p010)**2 + 4*(abs(bet)**2) ))/2
    lamb0m = (p100+p010)/2 - (np.sqrt( (p100-p010)**2 + 4*(abs(bet)**2) ))/2
    delt0 = (p100-p010)/2
    #if (abs(delt)< 1E-10):
        #th = np.pi/2
    #else:    
    th1=arctanaux(abs(alp),delt1)
    fase1 = cmath.phase(alp)
    fase = np.exp(1j*fase1)
    th0=arctanaux(abs(bet),delt0)
    fase0 = cmath.phase(bet)
    fase00 = np.exp(1j*fase0)
    #if (abs(alp) > 0 ):
     #   fase += alp/abs(alp)
    alp1 = (np.cos(th1/2)**2)*logaux(lamb1p) + (np.sin(th1/2)**2)*logaux(lamb1m)
    bet1 = (np.sin(th1/2)**2)*logaux(lamb1p) + (np.cos(th1/2)**2)*logaux(lamb1m)   
    c1 =  np.cos(th1/2)*np.sin(th1/2)*(fase)*(logaux(lamb1p) - logaux(lamb1m) )
    
    alp0 = (np.cos(th0/2)**2)*logaux(lamb0p) + (np.sin(th0/2)**2)*logaux(lamb0m)
    bet0 = (np.sin(th0/2)**2)*logaux(lamb0p) + (np.cos(th0/2)**2)*logaux(lamb0m)   
    c0 =  np.cos(th0/2)*np.sin(th0/2)*(fase00)*(logaux(lamb0p) - logaux(lamb0m) )
    #print("c0")
    #print(c0)
    #print("c1")
    #print(c1)
    #print("alp0")
    #print(alp0)
    #print("bet0")
    #print(bet0)
    return alp1,bet1,c1,alp0,bet0,c0

def IcoheL(rho,E,U,Uf,gl,glu,mul,betal):
    p111,p110,p101,p100 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p011,p010,p001,p000 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    alp,bet = rho[2,4],rho[3,5]
    alp1,bet1,c1,alp0,bet0,c0 = Log2(rho)
    aux0 = -gl*(1-fermi(E,mul,betal))*(p100*alp0 + (bet*c0.conjugate()).real )
    aux1 = gl*fermi(E,mul,betal)*(p000*alp0)
    aux2 = gl*(1-fermi(E+Uf,mul,betal))*(p110*bet0)
    aux3 = - gl*(fermi(E+Uf,mul,betal))*(p010*bet0 + (bet*c0.conjugate()).real)
    aux4 = -glu*(1-fermi(E+U,mul,betal))*(p101*alp1 + (alp*c1.conjugate()).real)
    aux5 = glu*fermi(E+U,mul,betal)*(p001*alp1)

    aux6 = glu*(1-fermi(E+U+Uf,mul,betal))*(p111*bet1) 
    aux7 = -glu*fermi(E+U+Uf,mul,betal)*(p011*bet1 + (alp*c1.conjugate()).real)
    total = aux0 + aux1 + aux2 +aux3 + aux4 +aux5+aux6+aux7
    return total

def IcoheR(rho,E,U,Uf,gr,gru,mur,betar):
    p111,p110,p101,p100 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p011,p010,p001,p000 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    alp,bet = rho[2,4],rho[3,5]
    alp1,bet1,c1,alp0,bet0,c0 = Log2(rho)
    aux0 = -gr*(1-fermi(E,mur,betar))*(p010*bet0 + (bet*c0.conjugate()).real )
    aux1 = gr*fermi(E,mur,betar)*(p000*bet0)
    aux2 = gr*(1-fermi(E+Uf,mur,betar))*(p110*alp0)
    aux3 = - gr*(fermi(E+Uf,mur,betar))*(p100*alp0 + (bet*c0.conjugate()).real)
    aux4 = -gru*(1-fermi(E+U,mur,betar))*(p011*bet1 + (alp*c1.conjugate()).real)
    aux5 = gru*fermi(E+U,mur,betar)*p001*bet1
    aux6 = gru*(1-fermi(E+U+Uf,mur,betar))*(p111*alp1) 
    aux7 = -gru*fermi(E+U+Uf,mur,betar)*(p101*alp1 + (alp*c1.conjugate()).real)
    total = aux0 + aux1 + aux2 +aux3 + aux4 +aux5+aux6+aux7
    return total

def IclassicL(rho,E,U,Uf,gl,glu,mul,betal):
    p111,p110,p101,p100 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p011,p010,p001,p000 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    bet,alp = rho[2,4],rho[3,5]
    aux0 = gl*(1-fermi(E,mul,betal))*(p100*logaux(p000))
    aux1 = -gl*fermi(E,mul,betal)*(p000*logaux(p000))
    aux2 = -gl*(1-fermi(E+Uf,mul,betal))*(p101*logaux(p110))
    aux3 =  gl*(fermi(E+Uf,mul,betal))*(p010*logaux(p110))
    aux4 = glu*(1-fermi(E+U,mul,betal))*(p101*logaux(p001))
    aux5 = -glu*fermi(E+U,mul,betal)*(p001*logaux(p001))#*bet1
    aux6 = -glu*(1-fermi(E+U+Uf,mul,betal))*(p111*logaux(p111))
    aux7 = glu*fermi(E+U+Uf,mul,betal)*(p011*logaux(p111))
    total = aux0 + aux1 + aux2 +aux3 + aux4 +aux5+aux6+aux7
    return total

def IclassicR(rho,E,U,Uf,gr,gru,mur,betar):
    p111,p110,p101,p100 = rho[0,0],rho[1,1],rho[2,2],rho[3,3]
    p011,p010,p001,p000 = rho[4,4],rho[5,5], rho[6,6], rho[7,7]
    bet,alp = rho[2,4],rho[3,5]
    #print(p000)
    aux0 = gr*(1-fermi(E,mur,betar))*(p010*logaux(p000))
    aux1 = -gr*fermi(E,mur,betar)*(p000*logaux(p000))
    aux2 = -gr*(1-fermi(E+Uf,mur,betar))*(p110*logaux(p110))
    aux3 =  gr*(fermi(E+Uf,mur,betar))*(p100*logaux(p110))
    aux4 = gru*(1-fermi(E+U,mur,betar))*(p011*logaux(p001))
    aux5 = -gru*fermi(E+U,mur,betar)*(p001*logaux(p001))#*bet1
    aux6 = -gru*(1-fermi(E+U+Uf,mur,betar))*(p111*logaux(p111))
    aux7 = gru*fermi(E+U+Uf,mur,betar)*(p101*logaux(p111))
    total = aux0+aux1+aux2+aux3+aux4+aux5+aux6+aux7
    return total

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
Ed = 0.
U0 = 40.
Uf = 500
g0 = 5/1000

eV = 6.5
mul1 = eV/2
mur1 = -eV/2
mud1 = 2*eV

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
betar,betad,betal = 1/100,1/100,1/100

#gr,grU = (1/100)*(1/6), 1/100
#gl,glU = 1/100, (1/100)*(1/6)
#gd,gdU = 1/50,1/50 
gr,grU = (1/100), 1/100
gl,glU = 1/100, (1/100)
gd,gdU = 1/100,1/100 
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



times = np.linspace(0,500,1000)
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
    #cal2 = DistR(H,Ls,Lr,Ll,rho0,ts,-2)
    #Probnt2.append(cal2.real)

plt.plot(times,Probnt1,label = "3 ocupados")
plt.plot(times,Probnt2, label = "ocupado r y l")
plt.plot(times,Probnt3, label = "ocupado l y d")
plt.plot(times,Probnt4, label = "ocupado l")
plt.scatter(times,Probnt5, label = "ocupado r y d")
plt.plot(times,Probnt6, label = "ocupado r")
plt.plot(times,Probnt7, label = "ocupado d")
plt.plot(times,Probnt8, label = "vacio total")
plt.legend()
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

eVs = np.linspace(0,1000,1000)
Ql = []
Qr = []
Qd = []
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
for ev in eVs:
    mud0 = 2
    U00 = 40 #10
    #mud0 = 1-U00/2
    #Con Ed0 = mud0 -U00/2,E0=4 hay flujo de energia pero un orden menor al de
    #flujo de informacion
    #Ed0 = 1
    Ed0 = mud0 -U00/2
    Uf0 = 500 #50
    #Probar condicion (U00/E0)<<1,Strasberg
    E0 = 4
    Ls0 = Dissipator(E0,Ed0,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,gl,glU,gr,grU,gd,gdU)
    H0 = Hamiltonian(E0,Ed0,U00,Uf0,g0)
    superop0 = Liouvillian(H0,Ls0)
    Ll0 = Dl(E0,U00,Uf0,ev/2,betal,gl,glU)
    Lr0 = Dr(E0,U00,Uf0,-ev/2,betar,gr,grU)
    Ld0 = Dd(Ed0,U00,mud0,betad,gd,gdU)
    Htd0 =  Htd(E0,Ed0,U00,Uf0)
    rhof = Propagate(rho0,superop,40000)
    Ql0,Qr0,Qd0,Sl0,Sr0,Sd0,El0,Er0,Ed0,Wl0,Wr0,Wd0,cohe0,concu0,Nl0,Nr0,Nd0 = currents(Htd0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,superop0,rho0,40000)
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


plt.plot(eVs,Ql,label = r'$\dot{Q}_{L}$')
plt.plot(eVs,Qr, label = r'$\dot{Q}_{R}$') 
plt.plot(eVs,Qd,label = r'$\dot{Q}_{d}$')
plt.plot(eVs,Nls,label = r'$\dot{N}_{L}$')
plt.plot(eVs,Nrs, label = r'$\dot{N}_{R}$') 
plt.plot(eVs,Nds,label = r'$\dot{N}_{d}$')
#plt.xscale("log")
plt.xlabel(r'$eV$',fontsize = 20)
plt.ylabel("Heat current",fontsize = 20)
plt.legend()
plt.show()   

plt.plot(eVs,Nls,label = r'$\dot{N}_{L}$', color = 'b')
plt.plot(eVs,Nrs, label = r'$\dot{N}_{R}$', color = 'r') 
plt.plot(eVs,Nds,label = r'$\dot{N}_{d}$', color = 'k')
#plt.xscale("log")
plt.xlabel(r'$eV$',fontsize = 20)
#plt.ylim(-0.0018, 0.0018) 
#plt.legend(loc='upper left')  
#plt.ylabel("Particle current",fontsize = 20)
plt.legend()
plt.show()  

plt.plot(eVs,Qlr, label = r'$\dot{Q}_{rl}$', color = 'b')
plt.plot(eVs,Qd,label = r'$\dot{Q}_{d}$', color = 'r')
plt.legend()
#plt.xscale("log")
plt.show()


plt.plot(eVs,Qd,label = r'$\dot{Q}_{d}$', color = 'b')
plt.plot(eVs,Id, label = r'$\dot{I}_{d}$', color = 'r')
plt.plot(eVs,Sds,label = r'$\dot{\sigma}_{d}$', color = 'k')
#plt.xscale("log")
plt.xlabel(r'$eV$',fontsize = 20)
plt.ylabel("Heat current",fontsize = 20)
plt.legend()
plt.show()   

plt.plot(eVs,Sls,label = r'$\dot{\sigma}_{L}$', color = 'b')
plt.plot(eVs,Srs, label = r'$\dot{\sigma}_{R}$', color = 'r') 
plt.plot(eVs,Sds,label = r'$\dot{\sigma}_{d}$', color = 'k')
plt.xlabel(r'$eV$',fontsize = 20)
plt.ylabel("Entropy production",fontsize = 20)
plt.legend()
plt.show()   

plt.plot(eVs,Slr)
plt.ylabel(r'$\dot{\sigma}_{LR}$',fontsize = 20, color = 'r')
plt.xlabel(r'eV', fontsize = 20)
plt.show()


plt.plot(eVs,Isl, label = r'$\dot{I}_{LR}$', color = 'b')
#plt.plot(eVs,Coher, label = r'$\mathcal{I}_{cohel}$')
#plt.plot(eVs,Cohel, label = r'$\mathcal{I}_{coher}$')
#plt.plot(eVs,Classl, label = r'$\mathcal{I}_{classL}$')
#plt.plot(eVs,Classr, label = r'$\mathcal{I}_{classR}$')
#plt.plot(eVs,cohesum, label = r'$\mathcal{I}_{coheLR}$')
#plt.plot(eVs, sumtot, linestyle='--', color='blue')
plt.plot(eVs,Id, label = r'$\dot{I}_{d}$', color = 'r')
plt.xlabel(r'$eV$',fontsize = 20)
plt.xscale("log")
plt.legend()
plt.show()


plt.plot(eVs,Ers, label = r'$\dot{E}_{R}$', color = 'b')
plt.plot(eVs,Els, label = r'$\dot{E}_{L}$', color = 'r')
plt.plot(eVs,Eds, label = r'$\dot{E}_{d}$', color = 'k')
plt.xlabel(r'$eV$',fontsize = 20)
#plt.xscale("log")
plt.legend()
plt.show()

plt.plot(eVs,Erl, label = r'$\dot{E}_{rl}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.plot(eVs,Qlr, label = r'$\dot{Q}_{rl}$')
plt.legend()
#plt.xscale("log")
plt.show()


plt.plot(eVs,Erl, label = r'$\dot{E}_{rl}$', color = 'b')
plt.plot(eVs,Isl, label = r'$\dot{I}_{rl}$', color = 'r')
plt.plot(eVs,Flr, label = r'$\dot{\mathcal{F}}_{rl}$', color = 'k')
plt.plot(eVs,Tisl,label = r'$T\dot{I}_{rl}$', color = 'c')
plt.plot(eVs,Wt,label = r'$\dot{W}_{LR}$', color = 'm')
plt.plot(eVs,Qlr,marker = 'D',linestyle='-',label = r'$\dot{Q}_{LR}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.legend()
#plt.xscale("log")
plt.show()

plt.plot(eVs,Eds, label = r'$\dot{E}_{d}$', color = 'b')
plt.plot(eVs,Isl, label = r'$\dot{I}_{rl}$', color = 'r')
plt.plot(eVs,Fd, label = r'$\dot{\mathcal{F}}_{d}$', color = 'k')
plt.plot(eVs,Tid,label = r'$T\dot{I}_{d}$', color = 'c')
plt.plot(eVs,Wdf,label = r'$\dot{W}_{d}$', color = 'm')
plt.plot(eVs,Qdf,marker = 'D',linestyle='-',label = r'$\dot{Q}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.legend()
#plt.xscale("log")
plt.show()


plt.plot(eVs,cohes, label = r'$\mathcal{C}_{l_{1}}$', color = 'b')
#plt.plot(eVs,Isl, label = r'$\dot{I}_{rl}$')
#plt.plot(eVs,Coher, label = r'$\mathcal{I}_{cohel}$')
#plt.plot(eVs,Cohel, label = r'$\mathcal{I}_{coher}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.legend(fontsize=16) 
#plt.xscale("log")
plt.show()

#ojo aqui, bajo eV=200, los puntos L y R parecen estar siendo medidos
#mientras que al superar esa vara L empieza a medir 
plt.plot(eVs,Id, label = r'$I_{d}$', color = 'b')
plt.plot(eVs,Ilf, label = r'$\dot{I}_{l}$', color = 'r')
plt.plot(eVs,Irf, label = r'$\dot{I}_{r}$', color = 'k')
plt.xlabel(r'$eV$',fontsize = 20)
plt.legend()
#plt.xscale("log")
plt.show()

plt.plot(eVs,cohes,label = r'$\mathcal{C}_{l_{1}}$', color = 'b')
plt.plot(eVs,concv, label = r'$\mathcal{C}_{on}$', color = 'r')  
plt.xlabel(r'$eV$',fontsize = 20)   
plt.legend()
plt.show()
