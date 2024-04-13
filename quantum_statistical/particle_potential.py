import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.animation import FFMpegWriter

#global constants#########################################################################################
hbar = 1.0
m = 1.0
x = np.linspace(-250,250,1000)
dx = (np.max(x) - np.min(x))/(len(x)-1)
t = np.linspace(0,500,1000)
dt = (np.max(t) - np.min(t))/(len(t)-1)
#############################################################################################################
class quantum_particle():
    def __init__(self, sigma, mean, energy, x,V):
        self.sigma = sigma
        self.mean = mean
        self.k = np.sqrt(2.0*m*energy/hbar**2.0)
        self.psi = np.zeros(len(x), dtype=np.complex128)
        self.ham = np.zeros((len(x), len(x)), dtype =np.complex128)
        self.identity = np.identity(len(x), dtype = np.complex128)
        
        #the wave form at t=0
        self.psi = (1.0/np.sqrt(self.sigma*np.sqrt(2.0*np.pi)))*np.exp(-(x- self.mean)**2.0/(2.0*self.sigma)**2.0)*np.exp(1j*self.k*x)

        for i in range(0, len(self.ham)):
            self.ham[i,i] = 2.0/(dx**2.0) + V[i]
            if i != len(self.ham)-1:
                self.ham[i+1,i] = -1.0/dx**2.0
                self.ham[i,i+1] = -1.0/dx**2.0

        #update A and B
        self.A = self.identity + 1j*dt*self.ham/2.0
        self.B = self.identity - 1j*dt*self.ham/2.0

        #update C
        self.C = np.matmul(np.linalg.inv(self.A),self.B)

    def time_evolve(self):
        #just update the particle's wave function
        self.psi = np.matmul(self.C,self.psi)

###################################################################################################################

V = np.zeros(len(x), dtype=np.complex128)
V[600:] = 1.8

marius = quantum_particle(5.0,0.0,1.0,x,V)
#print(np.sum(np.abs(marius.psi)**2.0*dx)) check first normalization



#####################################################################################################

fig = plt.figure()

fig.set_size_inches(18.5, 10.5)
carla = plt.subplot(1,1,1)

metadata = dict(title='Movie', artist ='ayush_roy')
writer = FFMpegWriter(fps=10, metadata = metadata)

with writer.saving(fig, 'particle_potential.mp4',100):
    for i in range(0,len(t)):
        marius.time_evolve()
    
        if i%8==0: 
            carla.plot(x,30*(np.abs(marius.psi)**2.0))
            carla.axvline(50,0,8,c='black',linestyle='dashed')
            
            carla.fill_between(x,0,8,where=x>50,facecolor='yellow',alpha=0.4)
            carla.annotate('V>0',(150,3.5),fontsize=20)
            
            carla.fill_between(x,0,8,where=x<50,facecolor='red',alpha=0.4)
            carla.annotate('V=0',(-150,3.5),fontsize=20)
            carla.set_title('Gaussian Wave Packet in Multiple Potentials')
            #plt.xlabel('X Coordinate')
            carla.set_xticks([])
            carla.set_ylabel('Probability Density Scaled by 30')
            carla.set_xlim([-250,250])
            carla.set_ylim([0,5.5])
            writer.grab_frame()
            plt.cla()
 #########################################################################################

        

