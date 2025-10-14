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

dp = (1/np.sqrt(2))*(dl + dr)
dm = (1/np.sqrt(2))*(dr - dl)
dpdag = (1/np.sqrt(2))*(dldag + drdag)
dmdag = (1/np.sqrt(2))*(drdag - dldag)

nd = np.matmul(dddag,dd)
np0 = np.matmul(dpdag,dp)
nm = np.matmul(dmdag,dm)

def Liouvillian( H,Ls, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    superL = sum( [np.kron(L.conjugate(),L) - 1/2 * (np.kron( np.eye(d), L.conjugate().T.dot(L)) +
                                                     np.kron( L.T.dot(L.conjugate()),np.eye(d) ))
                                                      for L in Ls ] )        
    return superH + superL

#operadores del disipador dd
def Dd(E,U,mud,betad,gammad):
    d = len(np0)
    auxd1 = np.sqrt( fermi(E,mud,betad)*gammad )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nm)),dddag )
    auxd2 = np.sqrt( (1-fermi(E,mud,betad))*gammad )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nm)),dd)
    auxd3 = np.sqrt( fermi(E+U,mud,betad)*gammad )*np.matmul( np.matmul((np.eye(d)-np0) ,nm) + np.matmul((np.eye(d)-nm) ,np0) ,dddag )
    auxd4 = np.sqrt( (1-fermi(E+U,mud,betad))*gammad )*np.matmul(np.matmul((np.eye(d)-np0) ,nm) + np.matmul((np.eye(d)-nm) ,np0),dd)
    auxd5 = np.sqrt( fermi(E+ (2*U),mud,betad)*gammad )*np.matmul( np.matmul(np0 ,nm),dddag )
    auxd6 = np.sqrt( (1-fermi(E+(2*U),mud,betad))*gammad )*np.matmul(np.matmul(np0 ,nm),dd)

    return [auxd1,auxd2,auxd3,auxd4,auxd5,auxd6]

#operadores del disipador positivos
def Dp(E,g,U,Uf,mul,betal,gammal):
    d = len(np0)
    Ep = E+g
    gl = gammal/2
    auxl1 = np.sqrt( fermi(Ep,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dpdag )
    auxl2 = np.sqrt( (1-fermi(Ep,mul,betal))*gl )*np.matmul( np.matmul((np.eye(d)-nm),(np.eye(d)-nd)),dp)
    auxl3 = np.sqrt( fermi(Ep+U,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nm) ,nd),dpdag )
    auxl4 = np.sqrt( (1-fermi(Ep+U,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-nm) ,nd),dp)
    auxl5 = np.sqrt( fermi(Ep+Uf,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nd) ,nm),dpdag )
    auxl6 = np.sqrt( (1-fermi(Ep+Uf,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-nd) ,nm),dp)
    auxl7 = np.sqrt( fermi(Ep+U+Uf,mul,betal)*gl )*np.matmul( np.matmul(nm,nd),dpdag )
    auxl8 = np.sqrt( (1-fermi(Ep+U+Uf,mul,betal))*gl )*np.matmul(np.matmul(nm,nd),dp)

    return [auxl1,auxl2,auxl3,auxl4,auxl5,auxl6,auxl7,auxl8]

#operadores del disipador dpositivos
def Dm(E,g,U,Uf,mul,betal,gammal):
    d = len(np0)
    Em = E-g
    gl = gammal/2
    auxl1 = np.sqrt( fermi(Em,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dmdag )
    auxl2 = np.sqrt( (1-fermi(Em,mul,betal))*gl )*np.matmul( np.matmul((np.eye(d)-np0),(np.eye(d)-nd)),dm)
    auxl3 = np.sqrt( fermi(Em+U,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-np0) ,nd),dmdag )
    auxl4 = np.sqrt( (1-fermi(Em+U,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-np0) ,nd),dm)
    auxl5 = np.sqrt( fermi(Em+Uf,mul,betal)*gl )*np.matmul( np.matmul((np.eye(d)-nd) ,np0),dmdag )
    auxl6 = np.sqrt( (1-fermi(Em+Uf,mul,betal))*gl )*np.matmul(np.matmul((np.eye(d)-nd) ,np0),dm)
    auxl7 = np.sqrt( fermi(Em+U+Uf,mul,betal)*gl )*np.matmul( np.matmul(np0,nd),dmdag )
    auxl8 = np.sqrt( (1-fermi(Em+U+Uf,mul,betal))*gl )*np.matmul(np.matmul(np0,nd),dm)

    return [auxl1,auxl2,auxl3,auxl4,auxl5,auxl6,auxl7,auxl8]

def Dissipator(E,g,Ed,U,Uf,mul,mur,mud,betal,betar,betad,gammal,gammar,gammad):
    DR = Dp(E,g,U,Uf,mur,betar,gammar) + Dm(E,g,U,Uf,mur,betar,gammar)
    DL = Dp(E,g,U,Uf,mul,betal,gammal) + Dm(E,g,U,Uf,mul,betal,gammal)
    DD = Dd(Ed,U,mud,betad,gammad)

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
    Ep,Em = E+g, E-g
    a1 = Ep*np0 + Em*nm + Ed*nd
    a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(np0,nd) +  np.matmul(nm,nd) ) + Uf*np.matmul(np0,nm) 
    return a1+a2+a3

def currents(Hs,mul,mur,mud,Ll,Lr,Ld,superop,rho0,t):
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

    Sl = -np.trace( np.matmul(Dl,aux) )
    Sr = -np.trace( np.matmul(Dr,aux) )
    Sd = -np.trace( np.matmul(Dd,aux) )

    El = np.trace( np.matmul( Dl,Hs  ) )
    Er = np.trace( np.matmul( Dr,Hs  ) )
    Ed = np.trace( np.matmul(Dd,Hs))

    return Ql.real, Qr.real, Qd.real, Sl.real, Sr.real,Sd.real, El.real, Er.real, Ed.real


E = 0
Ed = 0.
U0 = 0.
Uf = 50
g0 = 30/100

eV = 6.5
mul1 = eV/2
mur1 = -eV/2
mud1 = 2*eV

betar,betad,betal = 1,1,1
gr,gd,gl = 1/100,1/300,1/100
#gr,gd,gl = 1/100,0,1/100

Ll = Dp(E,g0,U0,Uf,mul1,betal,gl) +  Dm(E,g0,U0,Uf,mul1,betal,gl)
Lr = Dp(E,g0,U0,Uf,mur1,betar,gr) +  Dm(E,g0,U0,Uf,mur1,betar,gr)
Ld = Dd(Ed,U0,mud1,betad,gd)


Ls = Dissipator(E,g0,Ed,U0,Uf,mul1,mur1,mud1,betal,betar,betad,gl,gr,gd)
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


times = np.linspace(0,3000,1000)
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

values =  Propagate(rho0,superop,3000)

plt.imshow(values.imag)
plt.colorbar()
plt.show()
    
Ufs = np.linspace(7,40,100)
Us = np.linspace(1,6,50)

coheUf = []
coheU = []


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

#for ts in times:
 #   Ql0,Qr0,Qd0 = currents(Htd0,mul1,mur1,mud1,Ll,Lr,Ld,superop,rho0,ts)
  #  Ql.append(Ql0)
   # Qr.append(Qr0)
    #Qd.append(Qd0)


##plt.plot(times,Ql,label = r'$\dot{Q}_{L}$')
#plt.plot(times,Qr, label = r'$\dot{Q}_{R}$') 
#plt.plot(times,Qd,label = r'$\dot{Q}_{d}$')
#plt.show()   

eVs = np.linspace(0,50,500)
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

for ev in eVs:
    mud0 = 2
    U00 = 3
    #mud0 = 1-U00/2
    Ed0 = 1
    Uf0 = 20

    Ls0 = Dissipator(E,g0,Ed0,U00,Uf0,ev/2,-ev/2,mud0,betal,betar,betad,gl,gr,gd)
    H0 = Hamiltonian(E,Ed0,U00,Uf0,g0)
    superop0 = Liouvillian(H0,Ls0)
    Ll0 = Dp(E,g0,U00,Uf0,ev/2,betal,gl) + Dm(E,g0,U00,Uf0,ev/2,betal,gl)
    Lr0 = Dp(E,g0,U00,Uf0,-ev/2,betar,gr) + Dm(E,g0,U00,Uf0,-ev/2,betar,gr)
    Ld0 = Dd(Ed0,U00,mud0,betad,gd)
    Hs0 = Hamiltonian(E,g0,Ed0,U00,Uf0)

    Ql0,Qr0,Qd0,Sl0,Sr0,Sd0,El0,Er0,Ed0 = currents(Hs0,ev/2,-ev/2,mud0,Ll0,Lr0,Ld0,superop0,rho0,30000)
    Ql.append(Ql0)
    Qr.append(Qr0)
    Qd.append(Qd0)
    sigmal = Sl0 - betal*Ql0
    sigmar = Sr0 - betar*Qr0
    Sls.append(Sl0 - betal*Ql0)
    Srs.append(Sr0 - betar*Qr0)
    Sds.append(Sd0 - betad*Qd0)
    Slr.append( sigmal + sigmar )
    Isl.append(-Sl0 - Sr0)
    Id.append(-Sd0)
    Els.append(El0)
    Ers.append(Er0)
    Eds.append(Ed0)
    Erl.append(El0 + Er0)


#plt.plot(eVs,Ql,label = r'$\dot{Q}_{L}$')
#plt.plot(eVs,Qr, label = r'$\dot{Q}_{R}$') 
plt.plot(eVs,Qd,label = r'$\dot{Q}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.ylabel("Heat current",fontsize = 20)
plt.legend()
plt.show()   


plt.plot(eVs,Sls,label = r'$\dot{\sigma}_{L}$')
plt.plot(eVs,Srs, label = r'$\dot{\sigma}_{R}$') 
plt.plot(eVs,Sds,label = r'$\dot{\sigma}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.ylabel("Entropy production",fontsize = 20)
plt.legend()
plt.show()   

plt.plot(eVs,Slr)
plt.ylabel(r'$\dot{\sigma}_{LR}$',fontsize = 20)
plt.xlabel(r'eV', fontsize = 20)
plt.show()


plt.plot(eVs,Isl, label = r'$\dot{I}_{LR}$')
plt.plot(eVs,Id, label = r'$\dot{I}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.xscale("log")
plt.legend()
plt.show()


plt.plot(eVs,Ers, label = r'$\dot{E}_{R}$')
plt.plot(eVs,Els, label = r'$\dot{E}_{L}$')
plt.plot(eVs,Eds, label = r'$\dot{E}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.xscale("log")
plt.legend()
plt.show()

plt.plot(eVs,Erl, label = r'$\dot{E}_{d}$')
plt.xlabel(r'$eV$',fontsize = 20)
plt.xscale("log")
plt.show()