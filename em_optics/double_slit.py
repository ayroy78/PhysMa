from tkinter import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import tkinter.messagebox as tmsg
import numpy as np

import matplotlib.cm as cm
#########################################################################################################
#the mathematical engine
def calc_A(x,y,k,L,d, phi,w,t):
    x1 = np.sqrt((x + d/2.0 )**2.0 + y**2.0 + L**2.0)
    x2 = np.sqrt((x - d/2.0 )**2.0 + y**2.0 + L**2.0)
    return np.cos(k*x1 - w*t + phi) + np.cos(k*x2 - w*t)

c = 3.0*10**8
#class for the GUI
#############################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x700")
        #self.maxsize(1100,700)
        #self.minsize(1100,700)
        #title of our gui
        self.title("Double Slit Interference")
        
    def show_plot(self, *args):
        lam = self.lam.get()*10**-9
        k = 2.0*np.pi/lam
        d = self.d.get()*10**-2
        L = self.L.get()*10**-2
        phi = self.phi.get()

        w = c*k
        t = 0.0 

        max_order = 10
        th_draw = np.arcsin(max_order*lam/d)
        edge_draw = L*np.tan(th_draw)
        N_screen = 100
        x_screen = np.linspace(-edge_draw, edge_draw, N_screen)
        y_screen = np.linspace(-edge_draw, edge_draw, N_screen)

        X_SCREEN,Y_SCREEN = np.meshgrid(x_screen, y_screen)
        
        Amplitude = calc_A(X_SCREEN, Y_SCREEN, k,L,d,phi,w,t)
        Intensity = Amplitude**2.0
        
        self.ax.cla()
        self.ax.pcolormesh(X_SCREEN, Y_SCREEN, Intensity, cmap = cm.gray)
        self.ax.set_xlabel('Screen (m)')
        self.ax.set_ylabel('Screen (m)')
        self.canvas.draw()
        

    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Double Slit Interference").pack()

    def set_experiment(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Set Experiment").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Wavelength (nm)").pack(side = LEFT, padx=10)
        self.lam = Scale(f2, bg ='red', fg='black',from_=400.0, to=700.0, orient=HORIZONTAL, resolution=10.0, length=200)
        self.lam.set(400.0)
        self.lam.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Slit Separation (cm)").pack(side = LEFT, padx=10)
        self.d = Scale(f2, bg ='red', fg='black',from_=0.1, to=1.0, orient=HORIZONTAL, resolution=0.1)
        self.d.set(0.5)
        self.d.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Screen Distance (cm)").pack(side = LEFT, padx=10)
        self.L = Scale(f2, bg ='red', fg='black',from_=10.0, to=100.0, orient=HORIZONTAL, resolution=1.0, length=200)
        self.L.set(50.0)
        self.L.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Phase").pack(side = LEFT, padx=10)
        self.phi = Scale(f2, bg ='red', fg='black',from_=0.0, to=2.0*np.pi, orient=HORIZONTAL, resolution=0.1, length=200)
        self.phi.set(0.0)
        self.phi.pack(side = LEFT, padx=10)

        self.lam.config(command = self.show_plot)
        self.L.config(command = self.show_plot)
        self.d.config(command = self.show_plot)
        self.phi.config(command = self.show_plot)
                
    def make_plotframe(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(8,4, forward=True)
        self.ax = self.fig.add_subplot(1,1,1)
        
        self.fig.subplots_adjust(bottom=0.2)

        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f3, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)
#################################################################################################

#instance of out GUI class
double_slit = GUI()

#methods that give life to the GUI
double_slit.make_title()
double_slit.set_experiment()
double_slit.make_plotframe()
double_slit.show_plot()

#the infinite GUI loop
double_slit.mainloop()