from tkinter import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import tkinter.messagebox as tmsg
import numpy as np

#mathematical equations
############################################################################################
calc_tht = lambda ni,nt,thi : np.arcsin(ni*np.sin(thi)/nt)

calc_r_val_p = lambda ni, nt, thi,tht : (nt*np.cos(thi) - ni*np.cos(tht))/(nt*np.cos(thi) + ni*np.cos(tht)) 

calc_t_val_p = lambda ni,nt,thi,tht : 2.0*ni*np.cos(thi)/(ni*np.cos(tht) + nt*np.cos(thi))

calc_r_val_s = lambda ni,nt,thi, tht: (ni*np.cos(thi) - nt*np.cos(tht) )/(nt*np.cos(tht) + ni*np.cos(thi)) 

calc_t_val_s = lambda ni,nt,thi,tht : 2.0*ni*np.cos(thi)/(ni*np.cos(thi) + nt*np.cos(tht))

##############################################################################################

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
        self.title("Fresnel Equations")
        
    def show_plot(self, *args):
            polarization = self.polarization.get()
            ni = self.ni.get()
            nt = self.nt.get()
            thi = self.thi.get()*np.pi/180.0 # convert to radians

            #check if transmission is allowed
            if (ni*np.sin(thi)/nt >=1.0):
                #not allowed, set r and t manually
                r = 1.0
                t = 0.0
               
            else:
                #allowed, use the equations to find tht,r and t
                tht = calc_tht(ni,nt,thi)
                if polarization == 'p':
                    r = calc_r_val_p(ni,nt,thi,tht)
                    t = calc_t_val_p(ni,nt,thi,tht)
                
                elif polarization == 's':
                    r = calc_r_val_s(ni,nt,thi,tht)
                    t = calc_t_val_s(ni,nt,thi,tht)
                    

            #incidence always has to be drawn
            i_p = -np.tan(np.pi/2.0 - thi)*(-1.0) + 0.0
            
            self.ax.cla()
            self.ax.axhline(y=0, linestyle='dashed', c='black')
            self.ax.axvline(x=0, linestyle='dashed', c='black')
            self.ax.plot([0.0,-1.0],[0.0,i_p],c='red', linewidth=5.0, label='i')
            #if transmission is not zero, draw the transmitted ray with the correct linewidth
            if t != 0:
                t_p = - np.tan(np.pi/2.0 - tht)*(1.0) + 0.0
                self.ax.plot([0.0,1.0],[0.0,t_p],c='green', linewidth=5.0*t, label='t')
            #if reflection is not zero, draw the reflected ray with the correct linewidth
            if r != 0:
                r_p = np.tan(np.pi/2.0 - thi)*(1.0) + 0.0
                self.ax.plot([0.0,1.0],[0.0,r_p],c='orange', linewidth=5.0*r, label='r')
            self.ax.set_xlim([-1,1])
            self.ax.set_ylim([-1,1])
            self.ax.legend(loc=3)

            # Medium 1
            self.ax.fill_between([-1,1], y1=0, y2=1, color='blue', alpha=ni - 1.0)

            # Medium 2
            self.ax.fill_between([-1,1], y1=0, y2=-1, color='blue', alpha=nt - 1.0)

            self.ax.xaxis.set_visible(False)
            self.ax.yaxis.set_visible(False)

            self.ax.annotate('n2',(0.5,-0.5))
            self.ax.annotate('n1',(-0.5,0.5))
            
            self.canvas.draw()

    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Fresnel Equations").pack()

    def set_incidence(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Set Incident Ray, Media and Polarization").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="n1").pack(side = LEFT, padx=10)
        self.ni = Scale(f2, bg ='red', fg='black',from_=1.0, to=2.0, orient=HORIZONTAL, resolution=0.1)
        self.ni.set(1.5)
        self.ni.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="n2").pack(side = LEFT, padx=10)
        self.nt = Scale(f2, bg ='red', fg='black',from_=1.0, to=2.0, orient=HORIZONTAL, resolution=0.1)
        self.nt.set(1.0)
        self.nt.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Incidence Angle").pack(side = LEFT, padx=10)
        self.thi = Scale(f2, bg ='red', fg='black',from_=0.1, to=90.0, orient=HORIZONTAL, resolution=1.0)
        self.thi.set(10.0)
        self.thi.pack(side = LEFT, padx=10)

        self.polarization = StringVar()
        self.polarization.set("Radio")
        self.polarization.set("p")

        p = Radiobutton(f2, bg = "lightblue", fg='black',text="p-polarization",  variable=self.polarization, value="p")
        p.pack(side=LEFT, padx=10)

        s = Radiobutton(f2, bg = "lightblue", fg='black',text="s-polarization",  variable=self.polarization, value="s")
        s.pack(side=LEFT,padx=10)

        self.ni.config(command = self.show_plot)
        self.nt.config(command = self.show_plot)
        self.thi.config(command = self.show_plot)
        p.config(command = self.show_plot)
        s.config(command = self.show_plot)
                
    def make_plotframe(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(8,4, forward=True)
        self.ax = self.fig.add_subplot(1,1,1)
        
        #self.fig.subplots_adjust(bottom=0.2)

        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f3, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)
#################################################################################################

#instance of out GUI class
fresnel = GUI()

#methods that give life to the GUI
fresnel.make_title()
fresnel.set_incidence()
fresnel.make_plotframe()
fresnel.show_plot()

#the infinite GUI loop
fresnel.mainloop()