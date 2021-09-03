# Dispersive Readout of Multilevel Systems

# by Tara Murphy 

#get_ipython().run_line_magic('matplotlib', 'inline')

#Import all relevant libraries

import numpy as np
import matplotlib as plt
from numpy import linalg as LA
import cmath
import matplotlib.pyplot as plt


class Dispersive_Readout:
    def __init__(self, g1, g2, k1, k2, gc, w, w0, delta, gamma, so, T):
        self.g1 = g1; self.g2 = g2
        self.k1 = k1; self.k2 = k2; self.k = k1 + k2
        self.gc = gc
        self.w = w; self.w0 = w0; self.dw = w - w0
        self.delta = delta
        self.gamma = gamma
        self. so = so
        self.T = T
     
     
       # self.a = self.alpha(B, self.g2, self.g2)
       # self.b = self.beta(B, self.g1, self.g2)
     
    def alpha(self, B):
        a = 0.5*(self.g1+self.g2)*B;
        return a
     
    def beta(self, B):
        b = 0.5*(self.g1-self,g2)*B;
        return b
    
    def prob(self, energies):
        T = self.T
        Z = 0;
        prob = [];
        k = 1
        for i in range(0,len(energies)):
            Z += np.exp(-energies[i]/(k*T))
            prob.append(np.exp(-energies[i]/(k*T)))
        for i in range(0, len(energies)):
            prob[i] = prob[i]/Z
        return prob
    
    def Sz(self, B, e):
        return np.array([[0, self.delta, 0], [self.delta, -e, self.so], [0, self.so, -self.alpha(B)]])
    
    def evaluate(self,B, e): 
        g1 = self.g1
        g2 = self.g2 
        T = self.T
        so = self.so
        gc = self.gc
        gamma = self.gamma 
        w0 = self.w0
        w = self.w
        k1 = self.k1; k2 = self.k2; k = k1 + k2;
        
           #Obtain the Hamiltonian and Z matrix    
        Sz1 = self.Sz(B, e);
        sigmaz = np.array([[1,0,0], [0,-1,0], [0,0,1]]);
        dw = w0 - w
        #Determine eigenvalues and Basis
        [eigvals, eigvec] = LA.eig(Sz1)
        eigvals2 = (eigvals)
           
           
        ind1 = [x for x in range(len(eigvals)) if eigvals[x] == eigvals2[0] ]
        ind2 = [x for x in range(len(eigvals)) if eigvals[x] == eigvals2[1]] 
        ind3 = [x for x in range(len(eigvals)) if eigvals[x] == eigvals2[2]] 
        eigvec2 = [eigvec[ind1[0]], eigvec[ind2[0]], eigvec[ind3[0]]]
           
           #Determine Probabilities 
        p = self.prob(eigvals2)
           
           
           #Determine Z in old basis
        Z =np.matmul( np.matmul(eigvec2, sigmaz), LA.inv(eigvec2))
           #print(Z)
           #Difference in Energies and Probabilties
        dE12 = eigvals2[0] - eigvals2[1];
        dE23 = eigvals2[1] - eigvals2[2];
        dp12 = p[0] - p[1];
        dp23 = p[1] - p[2];
           
           #Z Matrix elements
        Z12 = Z[0,1]; Z23 = Z[1,2];
        Z21 = Z[1,0]; Z32 = Z[2,1];
          
           #Calculating Chi
        x12 = dp12*( (abs(Z12)**2)/(w + dE12 + 1j*gamma*0.5) - (abs(Z21)**2)/(w - dE12 + 1j*gamma*0.5) );
        x23 = dp23*( (abs(Z23)**2)/(w + dE23 + 1j*gamma*0.5) - (abs(Z23)**2)/(w - dE23 + 1j*gamma*0.5) );
        chi = x12 + x23;
          
        x= 1+(1j*( (k1) )/(self.dw-(1j*k/2) + (gc**2)*chi))
           
           #Calculating phase and transmission coefficient 
        S11 = pow(abs(x), 2)
        phase = cmath.phase(x)
           
        return [S11, phase, p[0], p[1], p[2]]
    
    def run(self):
    #Number of points you would like to run for detuning and magnetic field
        N = 150
           
           #Minimum and Maximum Values of Magnetic Field and Detuning 
        [Bmin, Bmax] =[0, 5] 
        [emin, emax] =[-3, 3] 
           
        B1 = np.linspace(Bmin, Bmax, N )
        e1 = np.linspace(emin, emax, N)
           
        z = [];
        ph = []
           
           #Create Heatmap 
           
        dw = self.w0 - self.w
        for B in B1:
            for e in e1:
                [S11, phase, p0, p1, p2] =  self.evaluate(B,e)
                z.append(S11)
                ph.append(phase)
                   
        z1 = ((np.reshape(z, (len(B1), len(e1)))))
        ph1 = ((np.reshape(ph, (len(B1), len(e1)))))
                   
       #All plotting completed below
       
        plt.figure()
           
        ax = plt.gca()
        plt.pcolormesh(e1, B1, z1, shading='auto')
        plt.colorbar()
           
        plt.xlabel('Detuning', fontsize = 10, fontweight="bold")
        plt.ylabel('Magnetic Field', fontsize = 10, fontweight="bold")
           
           #title = "$\omega_0$ = "+ str(w01)
           #plt.title(title, fontweight="bold", fontsize = 14)
           
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
           
           #Uncomment below if you want to save file    
           #title = "w0 = "+ str(w01)+"3x3.png"
           #plt.savefig(title)
           
        plt.show()
