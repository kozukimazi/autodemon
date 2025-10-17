import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.linalg import expm
from scipy import integrate
from scipy.linalg import eigh, eig

######################################
###########funciones##################
######################################

def fermi(E,mu,beta):
    return 1/(np.exp((E-mu)*beta) + 1)

#######################################
###########JordanWigner################
#######################################

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
Nop = nl+nr+nd
#########################################
##############Hamiltonian################
#########################################

def Hamiltonian(El,Er,Ed,U,Uf,g):
    a1 = El*nl + Er*nr + Ed*nd
    a2 = g*( np.matmul(dldag,dr) + np.matmul(drdag,dl) )
    a3 = U* (np.matmul(nl,nd) +  np.matmul(nr,nd) ) + Uf*np.matmul(nl,nr) 
    return a1+a2+a3

def Hamiltonianthermo(El,Er,Ed,U,Uf,g):
    a1 = El*nl + Er*nr + Ed*nd    
    a3 = U* (np.matmul(nl,nd) +  np.matmul(nr,nd) ) + Uf*np.matmul(nl,nr) 
    return a1+a3
#####################################
#############thermalization##########
#####################################
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

def thermal(Hs,beta):
    thermic = expm(-beta*Hs)
    canonic = thermic/np.trace(thermic)

    
    p000 = np.matmul( v00.T.conj(), np.matmul( canonic,v00 ) )
    p100 = np.matmul( v100.T.conj(), np.matmul( canonic,v100 ) )
    p010 = np.matmul( v010.T.conj(), np.matmul( canonic,v010 ) )
    p001 = np.matmul( v001.T.conj(), np.matmul( canonic,v001 ) )
    p110 = np.matmul( v110.T.conj(), np.matmul( canonic,v110 ) )
    p101 = np.matmul( v101.T.conj(), np.matmul( canonic,v101 ) )
    p011 = np.matmul( v011.T.conj(), np.matmul( canonic,v011 ) )
    p111 = np.matmul( v111.T.conj(), np.matmul( canonic,v111 ) )

    return p000[0,0],p100[0,0],p010[0,0],p001[0,0], p110[0,0],p101[0,0],p011[0,0],p111[0,0]



####################################
######twoleveldiagonalization#######
####################################

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


#construimos los operadores de salto para redfield
def ladderup(el,er,g,Uf,U):
    dup,dupdag,dmin,dmindag = twolevel(el,er,g)
    Delta = (el-er)/2
    Prom =  (el+er)/2
    eup = Prom + np.sqrt(Delta**2 + g**2)
    emin = Prom - np.sqrt(Delta**2 + g**2)
    nup = np.matmul(dupdag,dup)
    nmin = np.matmul(dmindag,dmin)
    op1 = np.matmul(dup,np.matmul( (np.eye(8) - nd), (np.eye(8) - nmin) ))
    op2 = np.matmul(dup,np.matmul( (np.eye(8) - nd),  nmin ))
    op3 = np.matmul(dup,np.matmul( (np.eye(8) - nmin),  nd ))
    op4 = np.matmul(dup,np.matmul(  nd,  nmin ))
    e1 = (eup)
    e2 = (eup+Uf)
    e3 = (eup + U)
    e4 = (eup + U + Uf)
    op1d = np.matmul(dupdag,np.matmul( (np.eye(8) - nd), (np.eye(8) - nmin) ))
    op2d = np.matmul(dupdag,np.matmul( (np.eye(8) - nd),  nmin ))
    op3d = np.matmul(dupdag,np.matmul( (np.eye(8) - nmin),  nd ))
    op4d = np.matmul(dupdag,np.matmul(  nd,  nmin ))
    e1d = eup
    e2d = eup+Uf
    e3d = eup + U
    e4d = eup + U + Uf
    listd1 = [op1,op2,op3,op4]
    enerd1 = [e1,e2,e3,e4]
    listd2 = [op1d,op2d,op3d,op4d]
    enerd2 = [e1d,e2d,e3d,e4d]
    return listd1,enerd1,listd2,enerd2

def laddermin(el,er,g,Uf,U):
    dup,dupdag,dmin,dmindag = twolevel(el,er,g)
    Delta = (el-er)/2
    Prom =  (el+er)/2
    eup = Prom + np.sqrt(Delta**2 + g**2)
    emin = Prom - np.sqrt(Delta**2 + g**2)
    nup = np.matmul(dupdag,dup)
    nmin = np.matmul(dmindag,dmin)
    op1 = np.matmul(dmin,np.matmul( (np.eye(8) - nd), (np.eye(8) - nup) ))
    op2 = np.matmul(dmin,np.matmul( (np.eye(8) - nd),  nup ))
    op3 = np.matmul(dmin,np.matmul( (np.eye(8) - nup),  nd ))
    op4 = np.matmul(dmin,np.matmul(  nd,  nup ))
    e1 = emin
    e2 = (emin+Uf)
    e3 = (emin + U)
    e4 = (emin + U + Uf)
    op1d = np.matmul(dmindag,np.matmul( (np.eye(8) - nd), (np.eye(8) - nup) ))
    op2d = np.matmul(dmindag,np.matmul( (np.eye(8) - nd),  nup ))
    op3d = np.matmul(dmindag,np.matmul( (np.eye(8) - nup),  nd ))
    op4d = np.matmul(dmindag,np.matmul(  nd,  nup ))
    e1d = emin
    e2d = emin+Uf
    e3d = emin + U
    e4d = emin + U + Uf
    listd1 = [op1,op2,op3,op4]
    enerd1 = [e1,e2,e3,e4]
    listd2 = [op1d,op2d,op3d,op4d]
    enerd2 = [e1d,e2d,e3d,e4d]
    return listd1,enerd1,listd2,enerd2

##aca empezamos a construir los disipadores L y R en funcion de d+ y d-
def ladderl(el,er,g,Uf,U):
    Delta = (el-er)/2
    if (Delta**2 + g**2>1E-9):
        aux = Delta/np.sqrt( Delta**2 + g**2)
        theta = math.acos(aux)
    else:
        theta = math.acos(1) 

    listup1,e1up,listup2,e2up=  ladderup(el,er,g,Uf,U)
    listmin1,e1min,listmin2,e2min=  laddermin(el,er,g,Uf,U)
    al = np.cos(theta/2)
    bl = -np.sin(theta/2)
    #aqui se le asigna cada operador en la lista su respectiva frecuencia de bohr
    list1 = [al*listup1[0],al*listup1[1],al*listup1[2],al*listup1[3],bl*listmin1[0],bl*listmin1[1],bl*listmin1[2],bl*listmin1[3]]
    ener1 = [e1up[0],e1up[1],e1up[2],e1up[3],e1min[0],e1min[1],e1min[2],e1min[3]]
    list2 = [al*listup2[0],al*listup2[1],al*listup2[2],al*listup2[3],bl*listmin2[0],bl*listmin2[1],bl*listmin2[2],bl*listmin2[3]]
    ener2 = [e2up[0],e2up[1],e2up[2],e2up[3],e2min[0],e2min[1],e2min[2],e2min[3]]
    
    return list1,ener1,list2,ener2

def ladderr(el,er,g,Uf,U):
    Delta = (el-er)/2
    if (Delta**2 + g**2>1E-9):
        aux = Delta/np.sqrt( Delta**2 + g**2)
        theta = math.acos(aux)
    else:
        theta = math.acos(1) 

    listup1,e1up,listup2,e2up=  ladderup(el,er,g,Uf,U)
    listmin1,e1min,listmin2,e2min=  laddermin(el,er,g,Uf,U)
    al = np.sin(theta/2)
    bl = np.cos(theta/2)
    #aqui la lista de los operadores posibles
    list1 = [al*listup1[0],al*listup1[1],al*listup1[2],al*listup1[3],bl*listmin1[0],bl*listmin1[1],bl*listmin1[2],bl*listmin1[3]]
    ener1 = [e1up[0],e1up[1],e1up[2],e1up[3],e1min[0],e1min[1],e1min[2],e1min[3]]
    list2 = [al*listup2[0],al*listup2[1],al*listup2[2],al*listup2[3],bl*listmin2[0],bl*listmin2[1],bl*listmin2[2],bl*listmin2[3]]
    ener2 = [e2up[0],e2up[1],e2up[2],e2up[3],e2min[0],e2min[1],e2min[2],e2min[3]]
    
    return list1,ener1,list2,ener2

def ladderd(ed,el,er,g,U):
    dup,dupdag,dmin,dmindag = twolevel(el,er,g)
    nup = np.matmul(dupdag,dup)
    nmin = np.matmul(dmindag,dmin)
    op1 = np.matmul(dd,np.matmul((np.eye(8)-nup) ,(np.eye(8)-nmin) ))
    op2 = np.matmul(dd,np.matmul(nup ,(np.eye(8)-nmin) ) + np.matmul(nmin ,(np.eye(8)-nup) ))
    op3 = np.matmul(dd,np.matmul(nup,nmin))
    e1 = (ed)
    e2 = (ed+U)
    e3 = (ed+(2*U))
    op1d = np.matmul(dddag,np.matmul((np.eye(8)-nup) ,(np.eye(8)-nmin) ))
    op2d = np.matmul(dddag,np.matmul(nup ,(np.eye(8)-nmin) ) + np.matmul(nmin ,(np.eye(8)-nup) ))
    op3d = np.matmul(dddag,np.matmul(nup,nmin))
    e1d = ed
    e2d = ed+U
    e3d = ed+(2*U)
    list1=[op1,op2,op3]
    ener1 = [e1,e2,e3]
    list2 = [op1d,op2d,op3d]
    ener2 = [e1d,e2d,e3d]
    return list1,ener1,list2,ener2

##le entrego una lista de gammas que coinciden con las energias
def Dd(betad,mud,gd,ed,el,er,g,U):
    list1,ener1,list2,ener2 = ladderd(ed,el,er,g,U)
    Num = len(list1)
    dim = len(nl)
    superD1 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    superD2 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    for j in range(Num):
        for jprim in range(Num):
            consta = gd*(1-fermi(ener1[j],mud,betad)+ 1-fermi(ener1[jprim],mud,betad))
            P1 = (consta/2)*(np.kron(list1[j].conjugate(), list1[jprim]))
            P2 = gd*((1-fermi(ener1[jprim],mud,betad))/2)*(np.kron(np.eye(dim),np.matmul(list1[j].transpose().conjugate(), list1[jprim])))
            P3 = gd*((1-fermi(ener1[j],mud,betad))/2)*(np.kron(np.matmul(list1[jprim].transpose(), list1[j].conjugate()),np.eye(dim) ))
            superD1 += (P1-P2-P3)
    
    for j in range(Num):
        for jprim in range(Num):
            consta = gd*(fermi(ener2[j],mud,betad)+fermi(ener2[jprim],mud,betad))
            P1 = (consta/2)*(np.kron(list2[j].conjugate(), list2[jprim]))
            P2 = gd*((fermi(ener2[jprim],mud,betad))/2)*(np.kron(np.eye(dim),np.matmul(list2[j].transpose().conjugate(), list2[jprim])))
            P3 = gd*((fermi(ener2[j],mud,betad))/2)*(np.kron(np.matmul(list2[jprim].transpose(), list2[j].conjugate()),np.eye(dim) ))
            superD2 += (P1-P2-P3)
    return superD1+superD2
##le entrego una lista de gammas que coinciden con las energias 
def Dlgen(betal,mul,gl1,el,er,g,Uf,U):
    list1,ener1,list2,ener2 = ladderl(el,er,g,Uf,U)
    Num = len(list1)
    dim = len(nl)
    superL1 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    superL2 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    for j in range(Num):
        for jprim in range(Num):
            consta = gl1[j]*(1-fermi(ener1[j],mul,betal))+ gl1[jprim]*(1 - fermi(ener1[jprim],mul,betal))
            P1 = (consta/2)*(np.kron(list1[j].conjugate(), list1[jprim]))
            P2 = gl1[jprim]*((1-fermi(ener1[jprim],mul,betal))/2)*(np.kron(np.eye(dim),np.matmul(list1[j].transpose().conjugate(), list1[jprim])))
            P3 = gl1[j]*((1-fermi(ener1[j],mul,betal))/2)*(np.kron(np.matmul(list1[jprim].transpose(), list1[j].conjugate()),np.eye(dim) ))
            superL1 += (P1-P2-P3)
    
    for j in range(Num):
        for jprim in range(Num):
            consta = (gl1[j]*fermi(ener2[j],mul,betal)+gl1[jprim]*fermi(ener2[jprim],mul,betal))
            P1 = (consta/2)*(np.kron(list2[j].conjugate(), list2[jprim]))
            P2 = gl1[jprim]*((fermi(ener2[jprim],mul,betal))/2)*(np.kron(np.eye(dim),np.matmul(list2[j].transpose().conjugate(), list2[jprim])))
            P3 = gl1[j]*((fermi(ener2[j],mul,betal))/2)*(np.kron(np.matmul(list2[jprim].transpose(), list2[j].conjugate()),np.eye(dim) ))
            superL2 += (P1-P2-P3)
    return superL1+superL2


def Drgen(betar,mur,gr1,el,er,g,Uf,U):
    list1,ener1,list2,ener2 = ladderr(el,er,g,Uf,U)
    Num = len(list1)
    dim = len(nl)
    superR1 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    superR2 = np.zeros((dim**2,dim**2) , dtype = np.complex_)
    for j in range(Num):
        for jprim in range(Num):
            consta = gr1[j]*(1 - fermi(ener1[j],mur,betar)) + gr1[jprim]*(1 -fermi(ener1[jprim],mur,betar))
            P1 = (consta/2)*(np.kron(list1[j].conjugate(), list1[jprim]))
            P2 = gr1[jprim]*((1-fermi(ener1[jprim],mur,betar))/2)*(np.kron(np.eye(dim),np.matmul(list1[j].transpose().conjugate(), list1[jprim])))
            P3 = gr1[j]*((1-fermi(ener1[j],mur,betar))/2)*(np.kron(np.matmul(list1[jprim].transpose(), list1[j].conjugate()),np.eye(dim) ))
            superR1 += (P1-P2-P3)
    
    for j in range(Num):
        for jprim in range(Num):
            consta = (gr1[j]*fermi(ener2[j],mur,betar)+gr1[jprim]*fermi(ener2[jprim],mur,betar))
            P1 = (consta/2)*(np.kron(list2[j].conjugate(), list2[jprim]))
            P2 = gr1[jprim]*((fermi(ener2[jprim],mur,betar))/2)*(np.kron(np.eye(dim),np.matmul(list2[j].transpose().conjugate(), list2[jprim])))
            P3 = gr1[j]*((fermi(ener2[j],mur,betar))/2)*(np.kron(np.matmul(list2[jprim].transpose(), list2[j].conjugate()),np.eye(dim) ))
            superR2 += (P1-P2-P3)
    return superR1+superR2



def superH( H, hbar = 1):
    d = len(H)
    superH = -1j/hbar * (np.kron(np.eye(d), H ) - np.kron(H.T,  np.eye(d))   )
    
    return superH 

def supertotal(H,Ll,Lr,Ld):
    supH = superH(H)
    return supH + Ll+Lr+Ld

def Propagate(rho0,superop,t):
    d = len(rho0)
    propagator = expm (superop *t)
    vec_rho_t = propagator @ np.reshape(rho0,(d**2,1))
    return np.reshape( vec_rho_t, (d,d) )

def propcurrent(rho0,superop,Nop1):
    d = len(rho0)
    propagator = (superop)
    vec_rho_t = propagator @ np.reshape(rho0,(d**2,1))
    rhof =  np.reshape( vec_rho_t, (d,d) )
    total = np.trace(np.matmul(Nop1,rhof) )
    return total

#######################################
############ejecucion##################
#######################################
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

#betal,betar,betad = 1/100,1/100,5/10
betal,betar,betad =1/100,1/100,1/2
eV = 245
mul,mur,mud = eV/2,-eV/2,2
gd=1/50
#gd = 1/100
#gl,glu = 1/100,(1/100)*(1/6)
#gr,gru = (1/100)*(1/6),1/100
gl,glu = 1/100,(1/100)*(1/6)
gr,gru = (1/100)*(1/6),1/100
el,er =0.,0
#g = 5/1000
g = 600
U,Uf = 40,500
ed=mud-(U/2)
#esta lista entrega un kappa por cada operador de salto de dl (son 8)
#gl1 = [gl,gl,glu,glu,gl,gl,glu,glu]
#gr1 = [gr,gr,gru,gru,gr,gr,gru,gru]
gl1 = [gl,gl,glu,gl,gl,gl,glu,gl]
gr1 = [gr,gr,gru,gr,gr,gr,gru,gr]


Num = 80000
g0s = np.linspace(0,1,Num)
cohev = []
concuv = []
g0m = []
for g0 in g0s:  
    #aqui hay funciones que no deberian estar para optimizarlo
    ev = 100
    print(g0)
    mud0 = 2
    U00 = 40 #10
    #mud0 = 1-U00/2
    #Con Ed0 = mud0 -U00/2,E0=4 hay flujo de energia pero un orden menor al de
    #flujo de informacion
    #Ed0 = 1
    Ed0 = mud0-U00/2
    Uf0 = 500 #50
    E0l = E0r = E0 = 0
    Ld = Dd(betad,mud,gd,Ed0,E0l,E0r,g0,U00)
    mul = ev/2
    mur = -ev/2
    Ll = Dlgen(betal,ev/2,gl1,E0l,E0r,g0,Uf0,U00)
    Lr = Drgen(betar,-ev/2,gr1,E0l,E0r,g0,Uf0,U00)
    Ham = Hamiltonian(E0l,E0r,Ed0,U00,Uf0,g0)
    Htd = Hamiltonianthermo(E0l,E0r,Ed0,U00,Uf0,g0)
    Superop =  supertotal(Ham,Ll,Lr,Ld)
    cal1f = Propagate(rho0,Superop,40000)
    cohev.append(abs(cal1f[5,3]) + abs(cal1f[4,2]) )
    cohesum = abs(cal1f[5,3] + cal1f[4,2])
    PD = cal1f[0,0].real + cal1f[1,1].real 
    P0 = cal1f[7,7].real + cal1f[6,6].real 
    concurrencef = 2*cohesum - 2*np.sqrt(P0*PD)
    if (concurrencef > 0):
        concuv.append(concurrencef)
    else:
        concuv.append(0)
    g0m.append(g0/gl)


archivo = open("redcoup","w")
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
    archivo.write("\n")