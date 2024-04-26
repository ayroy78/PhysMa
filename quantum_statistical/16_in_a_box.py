from matplotlib import pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import tkinter.messagebox as tmsg

#############################################################################################
box_width = 100.0
box_height = 100.0
dt = 1/10.0

#declare 16 particles and arrange them in a box
x_start = box_width/4.0
y_start = box_height/4.0
x_stop = box_width*0.75
y_stop = box_height*0.75
dx = (x_stop - x_start)/(4-1)
dy = (y_stop -y_start)/(4-1)
xpos = x_start
ypos = y_start

##################################################################################################
class free_particle():
    def __init__(self, mass,xinit,yinit):
        self.m = mass
        self.x0 = xinit
        self.y0 = yinit
        self.x1 = xinit
        self.y1 = yinit
        self.x = xinit
        self.y = yinit
        self.vx = 0.0
        self.vy = 0.0
        self.ax = 0.0
        self.ay = 0.0
        
    def move(self):
        '''
        xnew = 2.0*self.x1 - self.x0 + (self.ax)*dt**2.0 
        ynew = 2.0*self.y1 - self.y0 + (self.ay)*dt**2.0

        self.vx = (xnew - self.x0)/(2.0*dt)
        self.vy = (ynew - self.y0)/(2.0*dt)
        
        #shifting values for future time step calculations....
        self.x0 = self.x1
        self.x1 = xnew
        self.y0 = self.y1
        self.y1 = ynew
        
        if self.x1<0 or self.x1>box_width:
            self.x1 = self.x0
            self.vx = -self.vx
            
        if self.y1<0 or self.y1>box_height:
            self.y1 = self.y0
            self.vy = - self.vy
        '''
        self.vx = self.vx + self.ax*dt
        self.vy = self.vy + self.ay*dt

        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt
        
        if self.x<0 or self.x>box_width:
            self.x = self.x - self.vx*dt
            self.vx = -self.vx
            
        if self.y<0 or self.y>box_height:
            self.y = self.y - self.vy*dt
            self.vy = - self.vy
###################################################################################################

free_particles = []
for i in range(0,4): #move in x direction
    xpos = x_start + i*dx
    ypos = y_start
    for j in range(0,4): #move in y direction
        ypos = y_start + j*dy
        free_particles.append(free_particle(1.0,xpos,ypos))
#################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x700")
        #self.maxsize(1100,700)
        #self.minsize(1100,700)
        #title of our gui
        self.title("16 Particles in a Box")
        self.state = 0
        self.speeds = np.zeros(16)
        self.time  = 0.0

    def show_evolution(self,*args):
        if self.state == 0:
            self.state = 1
        elif self.state == 1:
            self.state = 0

        while self.state == 1:
            #find the accelerations for each particle; do not update any position yet 
            g = self.g.get()
            for i in free_particles:
                #zero the accelerations
                i.ax = 0.0
                i.ay = 0.0
                for j in free_particles:
                    if i!= j:
                        delta_x = i.x - j.x
                        delta_y = i.y - j.y
                        r2 = delta_x**2.0 + delta_y**2.0
                        i.ax = i.ax + (1.0/i.m)*(g/r2**1.5)*delta_x
                        i.ay = i.ay + (1.0/i.m)*(g/r2**1.5)*delta_y
            self.time = self.time + dt        
            #move particles
            for i in free_particles:
                i.move()

            #Show real time motion
        
            self.ax.cla()
            self.ax.set_title('Real Time Motion, t:' + str(int(self.time)))
            for i in free_particles:
                self.ax.scatter(i.x,i.y, c='green')
            self.ax.set_xlim([0,box_width])
            self.ax.set_ylim([0,box_height])

        

            #Update speed Histogram
            self.ax2.cla()
            self.ax2.set_title('Speed Distribution')
            for i in range(0,len(self.speeds)):
                self.speeds[i] = np.sqrt(free_particles[i].vx**2.0 + free_particles[i].vy**2.0)

            self.ax2.hist(self.speeds)
            self.ax2.set_xlabel('Speeds')
            self.ax2.set_ylabel('Counts')
        
            self.canvas.draw()
            self.canvas.start_event_loop(0.1)

    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "16 Particles in a Box").pack()

    def set_interaction(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Fix Interaction Strength").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="g").pack(side = LEFT, padx=10)
        self.g = Scale(f2, bg ='red', fg='black',from_=10.00, to=100.00, orient=HORIZONTAL, resolution=0.01, length=400)
        self.g.set(50.0)
        self.g.pack(side = LEFT, padx=10)
        
        play_pause = Button(f2, fg='black', bg='red', text="Play/Pause", font = ("calibiri", 12, "bold"), command = self.show_evolution)
        play_pause.pack(side = LEFT, padx=10)
        
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

################################################################################################

#instance of the GUI class
ensemble = GUI()

#methods of the GUI
ensemble.make_title()
ensemble.set_interaction()
ensemble.make_plotframe()
ensemble.show_evolution()

#infinite loop of the GUI
ensemble.mainloop()


