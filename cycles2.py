import numpy as np
import matplotlib.pyplot as plt
import scipy


class rates:
    def _init_(self,a):
        self.a = a
    @classmethod
    def FermiU(cls,T,E,mu,U,n):
        
        aux = (E-mu + n*U)/T   
        aux1 = np.exp(aux) + 1
        return  1/aux1 
    

    @classmethod
    def rate(cls,TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob):
        

        #Fermi functions
        f1y = cls.FermiU(TD,ex,muD,U,1)
        f0y = cls.FermiU(TD,ex,muD,U,0)

        f1xL = cls.FermiU(TL,ey,muL,U,1) 
        f0xL = cls.FermiU(TL,ey,muL,U,0) 

        f1xR = cls.FermiU(TR,ey,muR,U,1) 
        f0xR = cls.FermiU(TR,ey,muR,U,0) 

        #Transitions for X dot
        W10y1 = gamma*f1y
        W01y1 = gamma*(1-f1y)

        W10y0 = gamma*f0y
        W01y0 = gamma*(1-f0y)

        #Transitions for lower Y dot

        W01x0L = gamma0L*(1-f0xL)
        W01x1L = gamma1L*(1-f1xL)

        W01x0R = gamma0R*(1-f0xR)
        W01x1R = gamma1R*(1-f1xR)

        W10x0L = gamma0L*(f0xL)
        W10x1L = gamma1L*(f1xL)

        W10x0R = gamma0R*(f0xR)
        W10x1R = gamma1R*(f1xR)
        
        #Total Transitiom in Y dot
        W01x0 = W01x0L + W01x0R 
        W01x1 = W01x1L + W01x1R 
        
        W10x0 = W10x0L + W10x0R
        W10x1 = W10x1L + W10x1R



        #p(0,0), p(1,0), p(0,1), p(1,1) (convencion)
        #necesitamos definir el sistema

        
        a11 = -( W10x0 + W10y0 )
        a12 = W01y0 
        a13 = W01x0

        a22 = -(W10x1 + W01y0)
        a21 = W10y0
        a24 = W01x1

        a33 = -(W01x0 + W10y1)
        a31 = W10x0
        a34 = W01y1

        a44 = -(W01x1 + W01y1)
        a42 = W10x1
        a43 = W10y1

        Mat    = np.array([[a11,a12,a13,0],
                           [a21,a22,0,a24],
                           [a31,0,a33,a34],
                           [0,a42,a43,a44]])
        def System(y,t):
            dpdt = np.dot(Mat,y)
            return dpdt

               

        t = np.linspace(0,1000,5000)
        solution = scipy.integrate.odeint(System,prob,t)
        return t,solution
    
    @classmethod
    def currentuse(cls,TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob):    
        t, stationary = cls.rate(TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob)
        #Fermi functions
        f1y = cls.FermiU(TD,ex,muD,U,1)
        f0y = cls.FermiU(TD,ex,muD,U,0)

        f1xL = cls.FermiU(TL,ey,muL,U,1) 
        f0xL = cls.FermiU(TL,ey,muL,U,0) 

        f1xR = cls.FermiU(TR,ey,muR,U,1) 
        f0xR = cls.FermiU(TR,ey,muR,U,0) 

        #Transitions for X dot
        W10y1 = gamma*f1y
        W01y1 = gamma*(1-f1y)

        W10y0 = gamma*f0y
        W01y0 = gamma*(1-f0y)

        #Transitions for lower Y dot

        W01x0L = gamma0L*(1-f0xL)
        W01x1L = gamma1L*(1-f1xL)

        W01x0R = gamma0R*(1-f0xR)
        W01x1R = gamma1R*(1-f1xR)

        W10x0L = gamma0L*(f0xL)
        W10x1L = gamma1L*(f1xL)

        W10x0R = gamma0R*(f0xR)
        W10x1R = gamma1R*(f1xR)
        
        #Total Transitiom in Y dot
        W01x0 = W01x0L + W01x0R 
        W01x1 = W01x1L + W01x1R 
        
        W10x0 = W10x0L + W10x0R
        W10x1 = W10x1L + W10x1R
    
        p00 = stationary[-1,0]
        p10 = stationary[-1,1] 
        p01 = stationary[-1,2] 
        p11 = stationary[-1,3] 

        #Lin notation
        Jc1x = W01y0*p10 - W10y0*p00
        Jc2x =  W01y1*p11 - W10y1*p01
        Jc1g = (W10x0R + W10x0L)*p00 - (W01x0R + W01x0L)*p01
        Jc1y = W01x0R*p01 - W10x0R*p00
        Jc2y = W01x1R*p11 - W10x1R*p10 
        #Now the marginal probabilities
        py1 = p01 + p11
        py0 = p00 + p10
        px0 = p00 + p01
        px1 = p10 + p11

        #Condicional probabilities
        condp11 = p11/py1
        condp00 = p00/py0
        condp10 = p10/py0
        condp01 = p01/py1

        F = -np.log( (condp11*condp00)/(condp10*condp01) )


        return Jc1x,Jc2x,Jc1g,Jc1y,Jc2y,F
    

    @classmethod
    def current(cls,TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob):
        t, stationary = cls.rate(TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob)
        #Fermi functions
        f1y = cls.FermiU(TD,ex,muD,U,1)
        f0y = cls.FermiU(TD,ex,muD,U,0)

        f1xL = cls.FermiU(TL,ey,muL,U,1) 
        f0xL = cls.FermiU(TL,ey,muL,U,0) 

        f1xR = cls.FermiU(TR,ey,muR,U,1) 
        f0xR = cls.FermiU(TR,ey,muR,U,0) 

        #Transitions for X dot
        W10y1 = gamma*f1y
        W01y1 = gamma*(1-f1y)

        W10y0 = gamma*f0y
        W01y0 = gamma*(1-f0y)

        #Transitions for lower Y dot

        W01x0L = gamma0L*(1-f0xL)
        W01x1L = gamma1L*(1-f1xL)

        W01x0R = gamma0R*(1-f0xR)
        W01x1R = gamma1R*(1-f1xR)

        W10x0L = gamma0L*(f0xL)
        W10x1L = gamma1L*(f1xL)

        W10x0R = gamma0R*(f0xR)
        W10x1R = gamma1R*(f1xR)
        
        #Total Transitiom in Y dot
        W01x0 = W01x0L + W01x0R 
        W01x1 = W01x1L + W01x1R 
        
        W10x0 = W10x0L + W10x0R
        W10x1 = W10x1L + W10x1R


        p00 = stationary[-1,0]
        p10 = stationary[-1,1] 
        p01 = stationary[-1,2] 
        p11 = stationary[-1,3]    
        #print(p00,p10,p01,p11)

        #Jcy0 = (W10x0*p10 - W01x0*p00)
        
        #this works well
        #Jcy1 =  (W01x1*p10 - W10x1*p11)
        
        Je = -(W01x0R*p01 - W10x0R*p00)  - (W01x1R*p11 - W10x1R*p10)
        return Je    
    
    @classmethod
    def entropy(cls,TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob):
        Jc1x,Jc2x,Jc1g,Jc1y,Jc2y,F = cls.currentuse(TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob)    

        Je = Jc1y + Jc2y
        Deltamu = muL-muR
        Si = -Je*(Deltamu/TL) + Jc1g*( (U/TD) - (U/TL) )
        Irate = Jc1g*F
        Siy = -Je*(Deltamu/TL) + Jc1g*( F -(U/TL) )
        return Je,Si,Irate,Siy
 
prob = [1,0,0,0]
#TL,TR,TD,muL,muR,muD,ex,ey,gamma,gamma1L,gamma0L,gamma1R,gamma0R,U,prob
TL, TR, TD = 100,100,2
muL, muR,muD = 1.1,0.9, 0.2
ex,ey = 1,1
gamma,gamma1L,gamma0L,gamma1R,gamma0R = 100,0.5,1.5,1.5,0.5
t,solution = rates.rate(TL,TR,TD,  muL,muR,muD,  ex,ey, gamma,gamma1L,gamma0L,gamma1R,gamma0R, 0.1,prob)

plt.plot(t,solution[:,0])
plt.plot(t,solution[:,1])
plt.plot(t,solution[:,2])
plt.plot(t,solution[:,3])

print(solution[-1,0]+solution[-1,1]+solution[-1,2]+solution[-1,3])
print(solution[-1,3])

plt.show()    


Num = 200
eVs = np.linspace(0,1000,Num)
Je = []
Jem = []
Si = []
Is = []
Ism = []
Siy = []
for ev in eVs:
    U = 40
    mud0 = 2
    gamma0 = 1/50
    gamma1L0,gamma0L0 = 1/600,1/100
    gamma1R0,gamma0R0= 1/100,1/600
    ey0 = 4
    eD = mud0 - U/2
    valJ,valSi,valI,valSiy = rates.entropy(TL,TR,TD,  ev/2,-ev/2,mud0, eD,ey0, gamma0,gamma1L0,gamma0L0,gamma1R0,gamma0R0, U,prob)
    Je.append(-valJ)
    Jem.append(valJ)
    Si.append(-valSi)
    Is.append(-valI)
    Ism.append(valI)
    Siy.append(-valSiy)



plt.plot(eVs,Je, label = r'$\mathcal{J}_{e}$')
plt.plot(eVs,Jem, label = r'$\mathcal{J}_{e}$')
#plt.plot(eVs,Si, label = r'$\mathcal{S}_{i}$')
#plt.plot(eVs,Is, label = r'$\dot{\mathcal{I}}$')
#plt.plot(eVs,Siy,linestyle = '--', label = r'$\mathcal{S}^{Y}_{i}$')
#plt.ylabel(r'$\mathcal{J}_{e}$',fontsize=14)
plt.xlabel(r'$eV$',fontsize=14)
plt.legend()
plt.show()


plt.plot(eVs,Is, label = r'$\dot{\mathcal{I}}$')
plt.plot(eVs,Ism, label = r'$\dot{\mathcal{I}}$')
plt.xlabel(r'$eV$',fontsize=14)
plt.legend()
plt.show()


archivo = open("stochastic","w")
decimal_places = 7
total_width = 8
format_str = f"{{:.{decimal_places}f}}" 
#format_str = f"{{:{total_width}.{decimal_places}f}}"
for i in range(Num):
    archivo.write( format_str.format(eVs[i])) #guarda el grado del nodo
    #archivo.write(str(xs[i])) 
    archivo.write(" ") 
    #archivo.write(str(ys[i]))
    archivo.write( format_str.format(Je[i]))
    archivo.write(" ") 
    #archivo.write(str(ys[i]))
    archivo.write( format_str.format(Is[i]))
    archivo.write("\n")

