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
        self.title("TOV Polytropic Star Model")

    def refresh_plots(self):
        #set the star
        star.set_initial_data(self.polytropic_index.get(), self.central_pressure.get())
        star.RK4()
        radius = star.radius
        mass = star.mass
        f = interpolate.interp1d(star.r, star.P,kind='cubic')
        new_r = np.linspace(0,radius,100)
        new_P = f(new_r)
        #Internal Structure of the Star
        phi = np.linspace(0,2.0*np.pi,1000)
        x = []
        y = []
        z= []
        for i in range(0,len(new_r)):
            for j in range(0, len(phi)):
                x.append(new_r[i]*np.cos(phi[j]))
                y.append(new_r[i]*np.sin(phi[j]))
                z.append(new_P[i])
        self.ax.cla()
        self.ax.set_title('Internal Structure of the Star')
        self.ax.scatter(x,y,c=z, s=0.1)

        

        #Mass-Radius Relationship
        self.ax2.cla()
        
        if self.polytropic_index.get() == 1:
            self.ax2.plot(n1_data[:,0],n1_data[:,1])
        elif self.polytropic_index.get() == 1.5:
            self.ax2.plot(n15_data[:,0],n15_data[:,1])
        elif self.polytropic_index.get() == 2.0:
            self.ax2.plot(n2_data[:,0],n2_data[:,1])
        elif self.polytropic_index.get() == 2.5:
            self.ax2.plot(n25_data[:,0],n25_data[:,1])
        elif self.polytropic_index.get() == 3.0:
            self.ax2.plot(n3_data[:,0],n3_data[:,1])

        self.ax2.axvline(x=radius, linestyle='dashed', c='red')
        self.ax2.axhline(y=mass, linestyle='dashed', c='red')
        self.ax2.set_title('M-R Relationship')
        self.ax2.set_xlabel('Radius')
        self.ax2.set_ylabel('Mass')
        
        self.canvas.draw()

    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "TOV Polytropic Star Model").pack()

    def make_star(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Fix Star").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Polytropic Index").pack(side = LEFT, padx=10)
        self.polytropic_index = Scale(f2, bg ='lightblue', fg='black',from_=1.0, to=3.0, orient=HORIZONTAL, resolution=0.5, length=100)
        self.polytropic_index.set(1.0)
        self.polytropic_index.pack(side = LEFT, padx=10)
        

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Central Pressure").pack(side = LEFT, padx=10)
        self.central_pressure = Scale(f2, bg ='lightblue', fg='black',from_=0.001, to=1.0, orient=HORIZONTAL, resolution=0.001, length=700)
        self.central_pressure.set(0.005)
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

#####################################################################################################
#class for our polytropic star
class polytropic_star():
    def __init__(self):
        self.n = 0.0
        self.T = 0.0
        self.k = 1.0
        #arrays to hold stellar information
        self.dr = 0.01
        self.m = []
        self.P = []
        self.r = []
        #this is the central pressure
        self.rhoc = 0.0
        #this will be the final mass and radius
        self.mass = 0.0
        self.radius = 0.0

    #function to calculate density given pressure
    def calc_rho(self, P):
        return (P/self.k)**(1.0/self.T) + P/(self.T-1.0)
        
    #the coupled odes to solve   
    def m_slope(self, r,rho):
        return 4.0*np.pi*r**2.0*rho
    
    def P_slope(self, r, m, P, rho):
        return -(rho*m/r**2.0)*(1.0 + P/rho)*(1.0 + 4.0*np.pi*P*r**3.0/m)*(1.0 - 2.0*m/r)**(-1.0)
    
    #set initial data and the second point by hand to avoid discontinuity    
    def set_initial_data(self,polytropic_index, Pc):
        #zero the arrays each time a calculation has to be performed
        self.m = []
        self.P = []
        self.r = []

        self.n = polytropic_index
        self.T = 1.0 + 1.0/polytropic_index

        #the first point
        self.m.append(0.0)
        self.P.append(Pc)
        self.r.append(0.0)

        #hack for the second point
        self.P.append(Pc)
        self.r.append(self.dr)
        self.rhoc = self.calc_rho(Pc)
        self.m.append(self.m[0] + self.m_slope(self.r[1], self.rhoc)*self.dr)
        
    #RK-4 method to solve for mass and pressure of the star
    def RK4(self):
        i = 1
        while self.P[i].imag == 0 and self.P[i].real >0:
            rho = self.calc_rho(self.P[i])
            mk1 = self.dr*self.m_slope(self.r[i],rho)
            Pk1 = self.dr*self.P_slope(self.r[i],self.m[i],self.P[i], rho)
            m_d = self.m[i] + mk1/2.0
            P_d = self.P[i] + Pk1/2.0
            r_d = self.r[i] + self.dr/2.0

            if P_d.imag != 0 or P_d.real<0:
                break
            rho = self.calc_rho(P_d)
            mk2 = self.dr*self.m_slope(r_d,rho)
            Pk2 = self.dr*self.P_slope(r_d,m_d,P_d, rho)
            m_d = self.m[i] + mk2/2.0
            P_d = self.P[i] + Pk2/2.0
            r_d = self.r[i] + self.dr/2.0
            
            if P_d.imag != 0 or P_d.real<0:
                break

            rho = self.calc_rho(P_d)
            mk3 = self.dr*self.m_slope(r_d,rho)
            Pk3 = self.dr*self.P_slope(r_d,m_d,P_d, rho)
            m_d = self.m[i] + mk3
            P_d = self.P[i] + Pk3
            r_d = self.r[i] + self.dr

            if P_d.imag != 0 or P_d.real<0:
                break

            rho = self.calc_rho(P_d)
            mk4 = self.dr*self.m_slope(r_d,rho)
            Pk4 = self.dr*self.P_slope(r_d,m_d,P_d, rho)

            self.m.append(self.m[i] + (1.0/6.0)*(mk1 + 2.0*mk2 + 2.0*mk3 + mk4))
            self.P.append(self.P[i] + (1.0/6.0)*(Pk1 + 2.0*Pk2 + 2.0*Pk3 + Pk4))
            self.r.append(self.r[i] + self.dr)

            if self.P[i+1].imag != 0 or self.P[i+1].real<0:
                break

            i = i + 1
        self.mass = self.m[-1]
        self.radius = self.r[-1]
    
######################################################################################################

'''
#debugging for one star
NS = polytropic_star()
NS.set_initial_data(3.0,100.0)
NS.RK4()
print("The mass of the star is: " + str(NS.mass))
print("The radius of the star is: " + str(NS.radius))
print("The central density of the star is: " + str(NS.rhoc))
'''




'''
#Pre-Computing the Mass-Radius relationship calculation
star = polytropic_star() 
n = 2.5 #1.0,1.5,2.0,2.5,3.0
central_pressures = np.concatenate((np.linspace(0.001,0.01,100),np.linspace(0.01,0.1,100),np.linspace(0.1,1.0,100)))
m_r = np.zeros((len(central_pressures),2))
#masses = np.zeros(len(central_pressures))
#radii = np.zeros(len(central_pressures))
for i in range(0, len(central_pressures)):
    star.set_initial_data(n, central_pressures[i])
    star.RK4()
    m_r[i,0] = star.radius
    m_r[i,1] = star.mass
np.savetxt('n_2_5.txt',m_r)'''


#read in pre-computed data
n1_data = np.loadtxt('n_1.txt')
n15_data = np.loadtxt('n_1_5.txt')
n2_data = np.loadtxt('n_2.txt')
n25_data = np.loadtxt('n_2_5.txt')
n3_data = np.loadtxt('n_3.txt')

#instance of the GUI class
polytropic_stars = GUI()
star = polytropic_star() 

#methods of the GUI
polytropic_stars.make_title()
polytropic_stars.make_star()
polytropic_stars.make_plotframe()
polytropic_stars.refresh_plots()

#infinite loop of the GUI
polytropic_stars.mainloop()
