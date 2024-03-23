from matplotlib import pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import tkinter.messagebox as tmsg
from scipy import interpolate

######################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x700")
        #self.maxsize(1100,700)
        #self.minsize(1100,700)
        #title of our gui
        self.title("White Dwarf Fermi Gas Model")

    def refresh_plots(self):
        #set the star
        star.set_initial_data(self.chi.get(), self.central_pressure.get())
        star.RK4()
        radius = star.radius
        mass = star.mass
        f = interpolate.interp1d(star.r, star.rho,kind='cubic')
        new_r = np.linspace(0,radius,100)
        new_rho = f(new_r)
        #Internal Structure of the Star
        phi = np.linspace(0,2.0*np.pi,1000)
        x = []
        y = []
        z= []
        for i in range(0,len(new_r)):
            for j in range(0, len(phi)):
                x.append(new_r[i]*np.cos(phi[j]))
                y.append(new_r[i]*np.sin(phi[j]))
                z.append(new_rho[i])
        self.ax.cla()
        self.ax.set_title('Internal Structure of the Star')
        self.ax.scatter(x,y,c=z, s=0.1)

        

        #Mass-Radius Relationship
        self.ax2.cla()
        
        if self.chi.get() == 0.46:
            self.ax2.plot(chi_46_data[:,0],chi_46_data[:,1])
        elif self.chi.get() == 0.47:
            self.ax2.plot(chi_47_data[:,0],chi_47_data[:,1])
        elif self.chi.get() == 0.48:
            self.ax2.plot(chi_46_data[:,0],chi_46_data[:,1])
        elif self.chi.get() == 0.49:
            self.ax2.plot(chi_49_data[:,0],chi_49_data[:,1])
        elif self.chi.get() == 0.50:
           self.ax2.plot(chi_50_data[:,0],chi_50_data[:,1])

        self.ax2.axvline(x=radius, linestyle='dashed', c='red')
        self.ax2.axhline(y=mass, linestyle='dashed', c='red')
        self.ax2.set_title('M-R Relationship')
        self.ax2.set_xlabel('Radius')
        self.ax2.set_ylabel('Mass')
        
        self.canvas.draw()

    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "White Dwarf Fermi Gas Model").pack()

    def make_star(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Fix Star").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Chi").pack(side = LEFT, padx=10)
        self.chi = Scale(f2, bg ='red', fg='black',from_=0.46, to=0.50, orient=HORIZONTAL, resolution=0.01, length=100)
        self.chi.set(0.46)
        self.chi.pack(side = LEFT, padx=10)
        

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Central Pressure").pack(side = LEFT, padx=10)
        self.central_pressure = Scale(f2, bg ='red', fg='black',from_=-2.0, to=5.0, orient=HORIZONTAL, resolution=0.1, length=700)
        self.central_pressure.set(0.0)
        self.central_pressure.pack(side = LEFT, padx=10)

        refresh_star = Button(f2, fg='black', bg='red', text="Refresh Star", font = ("calibiri", 12, "bold"), command = self.refresh_plots)
        refresh_star.pack(side = LEFT, padx=10)
        
    def make_plotframe(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(10,4, forward=True)
        self.ax = self.fig.add_subplot(1,2,1)
        self.ax2 = self.fig.add_subplot(1,2,2)
        
        self.fig.subplots_adjust(bottom=0.2)
        self.fig.subplots_adjust(wspace=0.35)

        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f3, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)

######################################################################################################
#mathematical constants in units of solar mass and earth radii
G = 66.7
mh = 1.67*10**-57
me = 9.11*10**-61
c = 300.0

#####################################################################################################
#class for our white dwarf star
class white_dwarf():
    def __init__(self):
        self.chi = 0.0
        #arrays to hold stellar information
        self.dr = 0.001
        self.m = []
        self.rho = []
        self.r = []
        #this is the central density
        self.rhoc = 0.0
        #this is the density at which electrons turn relativistic
        self.rho0 = 0.0
        #this will be the final mass and radius
        self.mass = 0.0
        self.radius = 0.0

    #intermediate function
    def calc_rho0(self):
        self.rho0 = (9.8363*10**-4)/self.chi
          
    #the coupled odes to solve   
    def m_slope(self, r,rho):
        return 4.0*np.pi*rho*r**2.0
    
    def rho_slope(self, r, m, rho):
        x = (rho/self.rho0)**(1.0/3.0)
        alpha = (x**2.0)/(3.0*np.sqrt(1.0 + x**2.0))

        return -(mh*G*m*rho)/(self.chi*me*alpha*(c*r)**2.0)
    
    #set initial data and the second point by hand to avoid discontinuity    
    def set_initial_data(self,chi, rhoc):
        #zero the arrays each time a calculation has to be performed
        self.m = []
        self.rho = []
        self.r = []

        self.chi = chi
        self.calc_rho0()
        self.rhoc = self.rho0*(10**rhoc)
        
        #the first point
        self.m.append(0.0)
        self.rho.append(self.rhoc)
        self.r.append(0.0)

        #hack for the second point
        self.rho.append(self.rhoc)
        self.r.append(self.dr)
        #self.m.append(0.0)
        self.m.append(self.m[0] + self.m_slope(self.r[1], self.rhoc)*self.dr)
        
    #RK-4 method to solve for mass and density of the star
    def RK4(self):
        i = 1
        while self.rho[i].imag == 0 and self.rho[i].real >0:
            mk1 = self.dr*self.m_slope(self.r[i],self.rho[i])
            rhok1 = self.dr*self.rho_slope(self.r[i],self.m[i],self.rho[i])
            m_d = self.m[i] + mk1/2.0
            rho_d = self.rho[i] + rhok1/2.0
            r_d = self.r[i] + self.dr/2.0

            if rho_d.imag != 0 or rho_d.real<0:
                break

            mk2 = self.dr*self.m_slope(r_d,rho_d)
            rhok2 = self.dr*self.rho_slope(r_d,m_d, rho_d)
            m_d = self.m[i] + mk2/2.0
            rho_d = self.rho[i] + rhok2/2.0
            r_d = self.r[i] + self.dr/2.0
            
            if rho_d.imag != 0 or rho_d.real<0:
                break

            mk3 = self.dr*self.m_slope(r_d,rho_d)
            rhok3 = self.dr*self.rho_slope(r_d,m_d,rho_d)
            m_d = self.m[i] + mk3
            rho_d = self.rho[i] + rhok3
            r_d = self.r[i] + self.dr

            if rho_d.imag != 0 or rho_d.real<0:
                break

            mk4 = self.dr*self.m_slope(r_d,rho_d)
            rhok4 = self.dr*self.rho_slope(r_d,m_d,rho_d)

            self.m.append(self.m[i] + (1.0/6.0)*(mk1 + 2.0*mk2 + 2.0*mk3 + mk4))
            self.rho.append(self.rho[i] + (1.0/6.0)*(rhok1 + 2.0*rhok2 + 2.0*rhok3 + rhok4))
            self.r.append(self.r[i] + self.dr)

            if self.rho[i+1].imag != 0 or self.rho[i+1].real<0:
                break

            i = i + 1
        self.mass = self.m[-1]
        self.radius = self.r[-1]
    
    #Euler Method to find mass and radius of the star
    def euler(self):
        i = 1
        while self.rho[i].imag == 0 and self.rho[i].real >0:
            self.m.append(self.m[i] + self.dr*self.m_slope(self.r[i],self.rho[i]))
            self.rho.append(self.rho[i] + self.dr*self.rho_slope(self.r[i],self.m[i],self.rho[i]))
            self.r.append(self.r[i] + self.dr)

            if self.rho[i+1].imag != 0 or self.rho[i+1].real<0:
                break

            i = i + 1
        self.mass = self.m[-1]
        self.radius = self.r[-1]


    
    
######################################################################################################

'''
#debugging for one star
trial = white_dwarf()
trial.set_initial_data(0.5,5.0)
trial.RK4()
#trial.euler()
print("The mass of the star is: " + str(trial.mass))
print("The radius of the star is: " + str(trial.radius))
print('The critical density of the star is: ' + str(trial.rho0))
print("The central density of the star is: " + str(trial.rhoc))
print('The size of the iteration was: '+ str(len(trial.m)))

'''



'''
#Pre-Computing the Mass-Radius relationship calculation
star = white_dwarf() 
chi = 0.50
central_pressures = np.linspace(-2,5,100)
m_r = np.zeros((len(central_pressures),2))
for i in range(0, len(central_pressures)):
    star.set_initial_data(chi, central_pressures[i])
    star.RK4()
    m_r[i,0] = star.radius
    m_r[i,1] = star.mass
np.savetxt('0_5_0.txt',m_r)
'''


#read in pre-computed data
chi_46_data = np.loadtxt('0_4_6.txt')
chi_47_data = np.loadtxt('0_4_7.txt')
chi_48_data = np.loadtxt('0_4_8.txt')
chi_49_data = np.loadtxt('0_4_9.txt')
chi_50_data = np.loadtxt('0_5_0.txt')

#instance of the GUI class
white_dwarfs = GUI()
star = white_dwarf() 

#methods of the GUI
white_dwarfs.make_title()
white_dwarfs.make_star()
white_dwarfs.make_plotframe()
white_dwarfs.refresh_plots()

#infinite loop of the GUI
white_dwarfs.mainloop()
