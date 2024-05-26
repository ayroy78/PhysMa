from tkinter import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
import tkinter.messagebox as tmsg
import numpy as np

#####################################################################################################
#the mathematical engine

v_slope_photon = lambda l,r,M:  (l**2.0)/(r**3.0)  - (3*M*l**2.0)/(r**4.0)
v_slope_massive = lambda l,r,M:  (l**2.0)/(r**3.0)  - M/r**2.0 - (3*M*l**2.0)/(r**4.0)
t_slope = lambda e,r,M : e*(1.0 -(2*M/r))**(-1.0)
phi_slope = lambda l,r : l/r**2.0

V_eff_photon = lambda l,r,M : (l**2.0)/(2.0*r**2.0) - (M*l**2.0)/(r**3.0)

V_eff_massive = lambda l,r,M : (l**2.0)/(2.0*r**2.0) - M/r - (M*l**2.0)/(r**3.0)

tau = np.linspace(0,3000,10001)
dtau = (np.max(tau)- np.min(tau))/(len(tau)-1)

########################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x700")
        #self.maxsize(1100,700)
        #self.minsize(1100,700)
        #title of our gui
        self.title("Schwarzschild Geodesic Simulator")
        self.r_event_max = 0
        
##############################################################################################################################################
    def show_plots_initial_conditions(self, *args):

        if self.particle.get() == 'photon':
            e = np.sqrt(2.0*self.H.get())
        elif self.particle.get() =='massive':
            e = np.sqrt(2.0*self.H.get() + 1)
        l = self.l.get()
        M = self.mass.get()

        self.r_event_max = 2.0*M

        self.r.config(from_=self.r_event_max + 0.1, to=self.r_event_max*100.0)
    
        self.ax.cla()

        test_r = np.linspace(self.r_event_max+0.1, self.r_event_max*100, 1000)
        if self.particle.get() == 'photon':
            test_V_eff = V_eff_photon(l,test_r,M)
            V_eff_val = V_eff_photon(l, self.r.get(), M)
        elif self.particle.get() == 'massive':
            test_V_eff = V_eff_massive(l, test_r, M)
            V_eff_val = V_eff_massive(l, self.r.get(), M)
        self.ax.plot(test_r, test_V_eff)
        if self.particle.get() == 'photon':
            self.ax.axhline(y=e**2.0/2.0, linestyle='dashed', c='red', label='H')
        elif self.particle.get() == 'massive':
            self.ax.axhline(y=(e**2.0 -1.0)/2.0, linestyle='dashed', c='red', label='H')
        if self.particle.get() == 'photon':
            self.H.config(from_ = 0.0, to_ = np.max(test_V_eff)+1.0)
        elif self.particle.get() == 'massive':
            self.H.config(from_ = -0.499, to_ = np.max(test_V_eff)+1.0)

        self.ax.axvline(x=self.r.get(), linestyle='dashed', c='blue')
        self.ax.axhline(y = V_eff_val, linestyle ='dashed', c='blue', label=r'$V_{eff}(r)$')
        self.ax.scatter(self.r.get(), V_eff_val, c='blue')
        self.ax.set_xlabel(r'$r$')
        self.ax.set_ylabel(r'$V_{eff}$')
        self.ax.set_title('Potential Diagram')
        self.ax.legend()
        

        self.canvas.draw()        
        

###################################################################################################################################
    def ray_trace(self):
        #this is only activated when the 'compute trajectory' button is pressed. so even though it depends on
        #variables that are changing, it doesn't have any components that need to be constantly updated

        t = 0.0
        phi = 0.0
            
        x = np.zeros(len(tau))
        y = np.zeros(len(tau))
        #frst set the values to be used for ray tracing
        if self.particle.get() == 'photon':
            e = np.sqrt(2.0*self.H.get())
        elif self.particle.get() =='massive':
            e = np.sqrt(2.0*self.H.get() + 1)
        l = self.l.get()
        M = self.mass.get()

        self.r_event_max = 2.0*M

        r = self.r.get()

        x[0] = r *np.cos(phi)
        y[0] = r *np.sin(phi)
        
        v = self.v.get()
        if self.particle.get() == 'photon':
            speed_func = v_slope_photon
        elif self.particle.get() =='massive':
            speed_func = v_slope_massive

        r_d = 0.0
        v_d = 0.0

        Ham = np.ones(len(tau))
        break_index = 0

        # Create a meshgrid of spherical coordinates
        phi_circle = np.linspace(0, 2.0*np.pi,100)

        # Convert spherical coordinates to Cartesian coordinates
        x_circle = self.r_event_max * np.cos(phi_circle)
        y_circle = self.r_event_max * np.sin(phi_circle)

        #then check if ray tracing is possible. if so, do the RK4
        if self.particle.get() == 'photon':
            V_eff_val = V_eff_photon(l, self.r.get(), M)
        elif self.particle.get() == 'massive':
            V_eff_val = V_eff_massive(l, self.r.get(), M)
        
        if self.particle.get() == 'photon':
            temp_Ham = e**2.0/2.0
        elif self.particle.get() == 'massive':
            temp_Ham = (e**2.0 -1.0)/2.0

        if (temp_Ham - V_eff_val)<0:
            tmsg.showinfo("Invalid Initial Condition","The particle is in the forbidden region. Put the particle in the allowed region.")
        else:
            #now set the speeds
            v_t = t_slope(e,self.r.get(),M)
            v_phi = phi_slope(l, self.r.get())
            if v == 'plus':
                v = np.sqrt(2.0*(temp_Ham-V_eff_val))
            elif v == 'minus':
                v = - np.sqrt(2.0*(temp_Ham-V_eff_val))
            
            Ham[0] = v**2.0/2.0 + V_eff_val
               
            #Calculate the ray
            for i in range(0,len(tau)-1):
                #start with (r[i],v[i]) and employ RK-4 algorithm to find the 4 slopes
                v_k1 = (dtau)*speed_func(l,r,M)
                r_k1 = (dtau)*v
                t_k1 = (dtau)*v_t
                phi_k1 = (dtau)*v_phi

                r_d = r + r_k1/2.0
                v_d = v + v_k1/2.0
                v_t = t_slope(e,r_d,M)
                v_phi = phi_slope(l,r_d)

                if (r_d <= self.r_event_max):
                    r = 0.0
                    break_index = i+1
                    x[i+1] = r *np.cos(phi)
                    y[i+1] = r *np.sin(phi)
                    Ham[i+1] = Ham[i]
                    break

                v_k2 = (dtau)*speed_func(l,r_d,M)
                r_k2 = (dtau)*(v_d)
                t_k2 = (dtau)*v_t
                phi_k2 = (dtau)*v_phi

                r_d = r + r_k2/2.0
                v_d = v + v_k2/2.0
                v_t = t_slope(e,r_d,M)
                v_phi = phi_slope(l,r_d)

                if (r_d <= self.r_event_max):
                    r = 0.0
                    break_index = i+1
                    x[i+1] = r *np.cos(phi)
                    y[i+1] = r *np.sin(phi)
                    Ham[i+1] = Ham[i]
                    break

                v_k3 = (dtau)*speed_func(l,r_d,M)
                r_k3 = (dtau)*(v_d)
                t_k3 = (dtau)*v_t
                phi_k3 = (dtau)*v_phi

                r_d = r + r_k3
                v_d = v + v_k3
                v_t = t_slope(e,r_d,M)
                v_phi = phi_slope(l,r_d)

                if (r_d <= self.r_event_max):
                    r = 0.0
                    break_index = i+1
                    x[i+1] = r *np.cos(phi)
                    y[i+1] = r *np.sin(phi)
                    Ham[i+1] = Ham[i]
                    break

                v_k4 = (dtau)*speed_func(l,r_d,M)
                r_k4 = (dtau)*(v_d)
                t_k4 = (dtau)*v_t
                phi_k4 = (dtau)*v_phi

                #then add the average slope increment
                if (r + (1.0/6.0)*(r_k1 + 2*r_k2 + 2*r_k3 + r_k4) <= self.r_event_max):
                    r = 0.0
                    break_index = i+1
                    x[i+1] = r *np.cos(phi)
                    y[i+1] = r *np.sin(phi)
                    Ham[i+1] = Ham[i]
                    break
                
                v = v + (1.0/6.0)*(v_k1 + 2*v_k2 + 2*v_k3 + v_k4)
                r = r + (1.0/6.0)*(r_k1 + 2*r_k2 + 2*r_k3 + r_k4)
                t = t + (1.0/6.0)*(t_k1 + 2*t_k2 + 2*t_k3 + t_k4)
                phi = phi + (1.0/6.0)*(phi_k1 + 2*phi_k2 + 2*phi_k3 + phi_k4)
                v_t = t_slope(e,r,M)
                v_phi = phi_slope(l,r)

                x[i+1] = r * np.cos(phi)
                y[i+1] = r * np.sin(phi)
                
                if self.particle.get() == 'photon':
                    V_eff_val = V_eff_photon(l,r,M)
                elif self.particle.get() == 'massive':
                    V_eff_val = V_eff_massive(l,r,M)
                Ham[i+1] =  v**2.0/2.0 + V_eff_val 

            for i in range(0, len(Ham)):
                if Ham[i] == 0:
                    Ham[i] = 10**-15

            #plotting the trajectory
            if break_index != 0:
                for i in range(0, break_index+1):
                    self.ax2.cla()
                    self.ax2.plot(x_circle, y_circle, color='white')
                    self.ax2.scatter(x[i],y[i], c='red')
                    self.ax2.plot(x[0:i+1],y[0:i+1], linestyle='dashed', c='red', label='H= '+f"{Ham[i]:.3f}" )
                    self.ax2.set_xlim([min(np.min(x)-1,-self.r_event_max),max(np.max(x)+1,self.r_event_max)])
                    self.ax2.set_ylim([min(np.min(y)-1,-self.r_event_max),max(np.max(y)+1,self.r_event_max)])
                    #self.ax2.set_axis_off()
                    self.ax2.set_title('Trajectory')
                    self.ax2.legend(loc=1)

                    self.canvas.draw()
                    self.canvas.start_event_loop(0.1)
                    

            else:
                for i in range(0, len(tau)):
                    if i%100 == 0:
                        self.ax2.cla()

                        self.ax2.plot(x_circle, y_circle, color='white')
                        self.ax2.scatter(x[i],y[i], c='red')
                        self.ax2.plot(x[0:i+1],y[0:i+1], linestyle='dashed', c='red', label='H= '+f"{Ham[i]:.3f}" )
                        self.ax2.set_xlim([min(np.min(x)-1,-self.r_event_max),max(np.max(x)+1,self.r_event_max)])
                        self.ax2.set_ylim([min(np.min(y)-1,-self.r_event_max),max(np.max(y)+1,self.r_event_max)])
                        #self.ax2.set_axis_off()
                        self.ax2.set_title('Trajectory')
                        self.ax2.legend(loc=1)

                        self.canvas.draw()
                        self.canvas.start_event_loop(0.1)
                        

#############################################################################################################################################
    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Schwarzschild Geodesic Simulator").pack()
######################################################################################################################################
    def make_blackhole_particle(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Set Black Hole and Particle").pack()

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Mass").pack(side = LEFT, padx=10)
        self.mass = Scale(f2, bg ='red', fg='black',from_=0.5, to=1, orient=HORIZONTAL, resolution=0.1)
        self.mass.set(0.5)
        self.mass.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="H").pack(side = LEFT, padx=10)
        self.H = Scale(f2, bg ='red', fg='black',from_=-1.0, to=3.0, orient=HORIZONTAL, resolution=0.01, length=400)
        self.H.set(-0.02)
        self.H.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="l").pack(side = LEFT, padx=10)
        self.l = Scale(f2, bg ='red', fg='black',from_=-10, to=10, orient=HORIZONTAL, resolution=0.1)
        self.l.set(2.4)
        self.l.pack(side = LEFT, padx=10)

        self.particle = StringVar()
        self.particle.set("Radio")
        self.particle.set("massive")

        photon = Radiobutton(f2, bg = "lightblue", fg='black',text="Photon",  variable=self.particle, value="photon")
        photon.pack(side=LEFT, padx=10)

        massive = Radiobutton(f2, bg = "lightblue", fg='black',text="Massive",  variable=self.particle, value="massive")
        massive.pack(side=LEFT,padx=10)

        self.mass.config(command = self.show_plots_initial_conditions)
        self.H.config(command = self.show_plots_initial_conditions)
        self.l.config(command = self.show_plots_initial_conditions)
        photon.config(command = self.show_plots_initial_conditions)
        massive.config(command = self.show_plots_initial_conditions)

###################################################################################################
    def make_initial_conditions_frame(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        Label(f3, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Initial Conditions").pack()
        
        Label(f3, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="r").pack(side = LEFT, padx=10)
        self.r = Scale(f3, bg ='red', fg='black',from_=1.1, to=100.1, orient=HORIZONTAL, resolution=0.1, length=400)
        self.r.set(10.0)
        self.r.pack(side = LEFT, padx=10)

        self.v = StringVar()
        self.v.set("Radio")
        self.v.set("plus")

        v_plus = Radiobutton(f3, bg = "lightblue", fg='black',text="Positive Radial Speed",  variable=self.v, value="plus")
        v_plus.pack(side=LEFT, padx=10)

        v_minus = Radiobutton(f3, bg = "lightblue", fg='black',text="Negative Radial Speed",  variable=self.v, value="minus")
        v_minus.pack(side=LEFT,padx=10)
        
        ray_trace = Button(f3, fg='black', bg='red', text="Show Trajectory", font = ("calibiri", 12, "bold"), command = self.ray_trace)
        ray_trace.pack(side = LEFT, padx=10)

        self.r.config(command = self.show_plots_initial_conditions)
        v_plus.config(command = self.show_plots_initial_conditions)
        v_minus.config(command = self.show_plots_initial_conditions)
        
        
#######################################################################################################################################################
    def make_plotframe(self):
        f4 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f4.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='lightblue',linewidth=7)
        self.fig.set_size_inches(8,3.2, forward=True)
        self.ax = self.fig.add_subplot(1,2,1)
        self.ax2 = self.fig.add_subplot(1,2,2)
        self.ax2.set_facecolor('black')
        
        self.fig.subplots_adjust(bottom=0.2)

        self.canvas = FigureCanvasTkAgg(self.fig, master=f4)  
        self.canvas.get_tk_widget().pack(pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f4, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X)

       


############################################################################################################################################

#The main body of the program, which just calls all the functions from before

#instance of out GUI class
geodesic_simulator = GUI()

#methods that give life to the GUI
geodesic_simulator.make_title()
geodesic_simulator.make_blackhole_particle()
geodesic_simulator.make_initial_conditions_frame()
geodesic_simulator.make_plotframe()
geodesic_simulator.show_plots_initial_conditions()
geodesic_simulator.ray_trace()


#the infinite GUI loop
geodesic_simulator.mainloop()