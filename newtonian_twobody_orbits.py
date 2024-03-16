from tkinter import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import numpy as np

#####################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x640")
        #self.maxsize(1100,640)
        #self.minsize(1100,640)
        #title of our gui
        self.title("Two Body Newtonian Orbits")
######################################################################################################
    def ray_trace(self):
        #get parameters and set initial conditions
        q = self.massratio.get()
        a = self.a.get()
        e = self.e.get()

        m2 = 1.0
        m1 = q*m2
        M = m1 + m2
        mu = m1*m2/M
        a1 = a/(1.0 + q)
        a2 = a*q/(1.0 + q)

        if e<1:
            th = np.linspace(0.0,3.0*2.0*np.pi,900)
        elif e==1:
            th = np.linspace(0.0,np.pi-0.1,900)
        else:
            th = np.linspace(0.0, np.arccos(-1.0/e)-0.1,900)
            print(th[0:10], th[-10:])
        x1 = np.zeros(len(th))
        x2 = np.zeros(len(th))
        x3 = np.zeros(len(th))
        y1 = np.zeros(len(th))
        y2 = np.zeros(len(th))
        y3 = np.zeros(len(th))

        #this is just to extract E and L, placing the system's constituents at their perihelions. Actual positions found later.
        if e<1:
            r1_0 = a1*(1.0 - e)
            r2_0 = -a2*(1.0 - e)
            v_rel = np.sqrt(M*((2.0/(np.abs(r1_0) + np.abs(r2_0)))- (1.0/a)))
        elif e==1:
            r1_0 = a1
            r2_0 = -a2
            v_rel = np.sqrt(M*2.0/(np.abs(r1_0) + np.abs(r2_0)))
        else:
            r1_0 = a1*(e - 1.0)
            r2_0 = -a2*(e - 1.0)
            v_rel = np.sqrt(M*((2.0/(np.abs(r1_0) + np.abs(r2_0)))+ (1.0/a)))
        
        E = (1.0/2.0)*mu*v_rel**2.0 - M*mu/(np.abs(r1_0) + np.abs(r2_0))
        L = mu*v_rel*(np.abs(r1_0) + np.abs(r2_0))

        if e<1:
            r_plot = np.linspace(a*(1-e)-0.1,a*(1+e)+0.1,500)
        elif e==1:
            r_plot = np.linspace(a-0.1,2.0*a/(1+e*np.cos(max(th)))+0.1,500)
        else:
            r_plot = np.linspace(a*(e-1.0)-0.1,a*(e**2.0 -1.0)/(1+e*np.cos(max(th)))+0.1,500) 
        V_eff = lambda r: L**2.0/(2.0*mu*r**2.0) - mu*M/r 
        
        if e<1:
            r1 = a1*(1.0 - e**2.0)/(1.0 + e*np.cos(th))
            r2 = a2*(1.0 - e**2.0)/(1.0 + e*np.cos(th))
            r3 = a*(1.0 - e**2.0)/(1.0 + e*np.cos(th))
        elif e == 1:
            r1 = a1*2.0/(1.0 + e*np.cos(th))
            r2 = a2*2.0/(1.0 + e*np.cos(th))
            r3 = a*2.0/(1.0 + e*np.cos(th))
        else:
            r1 = a1*(e**2.0 - 1.0)/(1.0 + e*np.cos(th))
            r2 = a2*(e**2.0 - 1.0)/(1.0 + e*np.cos(th))
            r3 = a*(e**2.0 - 1.0)/(1.0 + e*np.cos(th))
        x1 = r1*np.cos(th)
        x2 = -r2*np.cos(th)
        y1 = r1*np.sin(th)
        y2 = -r2*np.sin(th)
        x3 = -r3*np.cos(th)
        y3 = -r3*np.sin(th)

        for i in range(0, len(th)):
            if i%9 == 0:
                #the effective potential
                self.ax2.cla()
                self.ax2.plot(r_plot, V_eff(r_plot), c='black')
                self.ax2.axhline(y=E, c='green', linestyle='dashed',label='E')
                self.ax2.scatter(r3[i],V_eff(r3[i]), c='green', label=r'$V_{eff}(r_\mu)$')
                self.ax2.set_title('Energy Diagram')
                self.ax2.set_ylabel(r'$V_{eff}(r)$')
                self.ax2.set_xlabel('r')
                self.ax2.legend(loc=4)

                #the trajectory
                self.ax.cla()
                self.ax.scatter(x1[i],y1[i], c='red', label='m1')
                self.ax.scatter(x2[i],y2[i], c='orange', label='m2')
                self.ax.scatter(x3[i],y3[i], c='green', label=r'$\mu$')
                self.ax.scatter(0,0,c='yellow', label='M')
                self.ax.plot(x1[0:i+1],y1[0:i+1], linestyle='dashed', c='red')
                self.ax.plot(x2[0:i+1],y2[0:i+1], linestyle='dashed', c='orange')
                self.ax.plot(x3[0:i+1],y3[0:i+1], linestyle='dashed', c='green')
                if e<1:
                    self.ax.set_xlim([min(np.min(x1)-1,np.min(x2)-1, np.min(x3)-1),max(np.max(x1)+1,np.max(x2)+1,np.max(x3)+1)])
                    self.ax.set_ylim([min(np.min(y1)-1,np.min(y2)-1, np.min(y3)-1),max(np.max(y1)+1,np.max(y2)+1,np.max(y3)+1)])
                else:
                    self.ax.set_xlim([min(np.min(x1[0:i+1])-1,np.min(x2[0:i+1])-1, np.min(x3[0:i+1])-1),max(np.max(x1[0:i+1])+1,np.max(x2[0:i+1])+1,np.max(x3[0:i+1])+1)])
                    self.ax.set_ylim([min(np.min(y1[0:i+1])-1,np.min(y2[0:i+1])-1, np.min(y3[0:i+1])-1),max(np.max(y1[0:i+1])+1,np.max(y2[0:i+1])+1,np.max(y3[0:i+1])+1)])
                #self.ax2.set_axis_off()
                self.ax.axhline(y=0,c='yellow')
                self.ax.axvline(x=0, c='yellow')
                self.ax.set_title('Trajectories')
                self.ax.legend(loc=4)

                self.canvas.draw()
                self.canvas.start_event_loop(0.1)
 
#############################################################################################################################################
    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Two Body Newtonian Orbits").pack()
######################################################################################################################################
    def make_binary_system(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Set Binary System Properties").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Mass Ratio (q=m1/m2)").pack(side = LEFT, padx=10)
        self.massratio = Scale(f2, bg ='red', fg='black',font ="bold", from_=1, to=100, orient=HORIZONTAL, resolution=1)
        self.massratio.set(1.0)
        self.massratio.pack(side=LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Semimajor axis (a)").pack(side = LEFT, padx=10)
        self.a = Scale(f2, bg ='red', fg='black',font ="bold",from_=2, to=10, orient=HORIZONTAL, resolution=1)
        self.a.set(1)
        self.a.pack(side=LEFT, padx=10)
        
        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Eccentricity (e)").pack(side = LEFT, padx=10)
        self.e = Scale(f2, bg ='red', fg='black',font="bold",from_=0., to=1.5, orient=HORIZONTAL, resolution=0.1)
        self.e.set(0.5)
        self.e.pack(side=LEFT, padx=10)
        
        ray_trace = Button(f2, fg='black', bg='red', text="Show Trajectory", font = ("calibiri", 12, "bold"), command = self.ray_trace)
        ray_trace.pack(side = LEFT, padx=70)

#######################################################################################################################################################
    def make_plotframe(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(10,4, forward=True)
        self.ax = self.fig.add_subplot(1,2,2)
        self.ax2 = self.fig.add_subplot(1,2,1)
        self.ax.set_facecolor('black')
        self.fig.subplots_adjust(bottom=0.15)
        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f3, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)

############################################################################################################################################

#The main body of the program, which just calls all the functions from before

#instance of out GUI class
binary = GUI()

#methods that give life to the GUI
binary.make_title()
binary.make_binary_system()
binary.make_plotframe()
binary.ray_trace()

#the infinite GUI loop
binary.mainloop()