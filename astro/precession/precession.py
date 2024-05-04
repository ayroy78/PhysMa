import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter
#####################################################################################
# Load the images
bg = mpimg.imread('bg.jpg')
aries = mpimg.imread('aries.jpg')
taurus = mpimg.imread('taurus.jpg')
gemini = mpimg.imread('gemini.jpg')
cancer = mpimg.imread('cancer.jpg')
leo = mpimg.imread('leo.jpg')
virgo = mpimg.imread('virgo.jpg')
libra = mpimg.imread('libra.jpg')
scorpio = mpimg.imread('scorpio.jpg')
sagittarius = mpimg.imread('sagittarius.jpg')
capricorn = mpimg.imread('capricorn.jpg')
aquarius = mpimg.imread('aquarius.jpg')
pisces = mpimg.imread('pisces.jpg')

######################################################################################
# Animation setup
fig, ax = plt.subplots()
R_zodiac = 100
phi = np.linspace(0, 2.0*np.pi,100)
x_zodiac = R_zodiac*np.cos(phi)
y_zodiac = R_zodiac*np.sin(phi)

R_earth = 50
x_earth = R_earth*np.cos(phi)
y_earth = R_earth*np.sin(phi)

rho = 30

phi0 = -np.pi/3.0
t0 = 0.0
frames = 400
omega = 2.0*np.pi/frames
speed = 25772.0/frames
############################################################################################
metadata = dict(title='Movie', artist ='ayush_roy')
writer = PillowWriter(fps=10, metadata = metadata)

with writer.saving(fig, 'precession.gif',100):
    for k in range(frames+1):

        # Plot each image at specified coordinates
        ax.imshow(bg, extent =[-105,105,-105,105])
        ax.imshow(aries, extent=[-30, -5, 70, 95])  # (left, right, bottom, top)
        ax.imshow(taurus, extent=[-60, -35, 40, 65])
        ax.imshow(gemini, extent=[-90, -65, 5, 35])
        ax.imshow(cancer, extent=[-90, -65, -35, -5])  # (left, right, bottom, top)
        ax.imshow(leo, extent=[-60, -35, -65, -40])
        ax.imshow(virgo, extent=[-30, -5, -95, -70])
        ax.imshow(libra, extent=[5, 30, -95, -70])  # (left, right, bottom, top)
        ax.imshow(scorpio, extent=[35, 60, -65, -40])
        ax.imshow(sagittarius, extent=[65, 90, -35, -5])
        ax.imshow(capricorn, extent=[65, 90, 5, 30])  # (left, right, bottom, top)
        ax.imshow(aquarius, extent=[40, 65, 40, 65])
        ax.imshow(pisces, extent=[5, 35, 65, 95])

        #The dividing white lines
        ax.plot(x_zodiac,y_zodiac, c='white', linestyle='dashed')
        ax.plot(x_earth,y_earth, c='white', linestyle='dashed')
        ax.plot([-100,-50],[0,0], c='white', linestyle='dashed')
        ax.plot([-R_zodiac*np.cos(np.pi/6),-R_earth*np.cos(np.pi/6)],[R_zodiac*np.sin(np.pi/6),R_earth*np.sin(np.pi/6)], c='white', linestyle='dashed')
        ax.plot([-R_zodiac*np.cos(np.pi/3),-R_earth*np.cos(np.pi/3)],[R_zodiac*np.sin(np.pi/3),R_earth*np.sin(np.pi/3)], c='white', linestyle='dashed')
        ax.plot([0,0],[50,100],c='white', linestyle='dashed')
        ax.plot([R_zodiac*np.cos(np.pi/3),R_earth*np.cos(np.pi/3)],[R_zodiac*np.sin(np.pi/3),R_earth*np.sin(np.pi/3)], c='white', linestyle='dashed')
        ax.plot([R_zodiac*np.cos(np.pi/6),R_earth*np.cos(np.pi/6)],[R_zodiac*np.sin(np.pi/6),R_earth*np.sin(np.pi/6)], c='white', linestyle='dashed')
        ax.plot([50,100],[0,0], c='white', linestyle='dashed')
        ax.plot([R_zodiac*np.cos(np.pi/6),R_earth*np.cos(np.pi/6)],[-R_zodiac*np.sin(np.pi/6),-R_earth*np.sin(np.pi/6)], c='white', linestyle='dashed')
        ax.plot([R_zodiac*np.cos(np.pi/3),R_earth*np.cos(np.pi/3)],[-R_zodiac*np.sin(np.pi/3),-R_earth*np.sin(np.pi/3)], c='white', linestyle='dashed')
        ax.plot([0,0],[-100,-50],c='white', linestyle='dashed')
        ax.plot([-R_zodiac*np.cos(np.pi/6),-R_earth*np.cos(np.pi/6)],[-R_zodiac*np.sin(np.pi/6),-R_earth*np.sin(np.pi/6)], c='white', linestyle='dashed')
        ax.plot([-R_zodiac*np.cos(np.pi/3),-R_earth*np.cos(np.pi/3)],[-R_zodiac*np.sin(np.pi/3),-R_earth*np.sin(np.pi/3)], c='white', linestyle='dashed')
        
        #The sun
        ax.scatter(0,0,c='yellow', s=800)

        #Update the earth's polar angle and the real physical time
        phi_earth = phi0 - omega*k
        time = t0 + speed*k

        #The earth
        ax.scatter(rho*np.cos(phi_earth), rho*np.sin(phi_earth), c='blue', s=200)

        #The line from earth to the zodiac
        ax.plot([rho*np.cos(phi_earth), -R_earth*np.cos(phi_earth)],[rho*np.sin(phi_earth), -R_earth*np.sin(phi_earth)], c='red', linestyle='dashed')

        #The time
        ax.annotate(str(int(time))+ ' y',(65,-95),c='white')

        ax.set_title('Precession: Shift of the Spring Equinox')
        ax.set_xlim([-105,105])
        ax.set_ylim([-105,105])
        ax.set_xticks([])
        ax.set_yticks([])
        
        writer.grab_frame()
        plt.cla()

###################################################################################################