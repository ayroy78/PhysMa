import pygame
import numpy as np

# Initialize Pygame
pygame.init()
screen_width, screen_height = 800, 600
screen = pygame.display.set_mode((screen_width, screen_height))
clock = pygame.time.Clock()
fps = 60

# Colors
white = (255, 255, 255)
black = (0, 0, 0)
red = (255,0,0)
green = (0,255,0)
blue = (0,0,255)

#The Pendulum Class################################################################################
g = 9.8
L = screen_height/4.0
w = np.sqrt(g/L)
r = 20.0
m = 1.0

class pendulum():
    def __init__(self,w,r,m):
        #natural frequency which can never change
        self.w = w
        #the radius of the pendulum bob
        self.radius = r
        #the mass of the pendulum bob
        self.mass = m
        #initial parameters
        self.th = np.pi/4.0
        self.t = 0.0
        self.v = 0.0
        self.dt = (1.0/60.0)
        self.x = L*np.sin(self.th)
        self.y = L*np.cos(self.th)
        self.E = (1.0/2.0)*self.mass*(self.v*L)**2.0 - self.mass*g*self.y
        #score
        self.score = 0

    def a(self,v,th,t,b,F,wext):
        return -b*v - self.w**2.0*np.sin(th) + F*np.cos(wext*t)
    
    def RK4(self,b,F,wext):
        vk1 = self.dt*self.a(self.v,self.th,self.t,b,F,wext)
        thk1 = self.dt*self.v

        vd = self.v + vk1/2.0
        thd = self.th + thk1/2.0
        td = self.t + self.dt/2.0

        vk2 = self.dt*self.a(vd,thd,td,b,F,wext)
        thk2 = self.dt*vd

        vd = self.v + vk2/2.0
        thd = self.th + thk2/2.0
        td = self.t + self.dt/2.0

        vk3 = self.dt*self.a(vd,thd,td,b,F,wext)
        thk3 = self.dt*vd

        vd = self.v + vk3
        thd = self.th + thk3
        td = self.t + self.dt

        vk4 = self.dt*self.a(vd,thd,td,b,F,wext)
        thk4 = self.dt*vd

        self.v = self.v + (1.0/6.0)*(vk1 + 2.0*vk2 + 2.0*vk3 + vk4)
        self.th = self.th + (1.0/6.0)*(thk1 + 2.0*thk2 + 2.0*thk3 + thk4)
        self.t = self.t + self.dt

        self.x = L*np.sin(self.th)
        self.y = L*np.cos(self.th)
        self.E = (1.0/2.0)*self.mass*(L*self.v)**2.0 - self.mass*g*self.y
       

    def draw(self):
        pygame.draw.circle(screen, red, (screen_width/2 + self.x, screen_height/2 + self.y),self.radius)
        

####################################################################################################################################

#The Bullet Class#############################################################################################################################
class bullet():
    def __init__(self,x,y,vx):
        #initial x and y position on the screen + speed
        self.x = x
        self.y = y
        self.vx = vx
    
    def move(self):
        self.x = self.x + self.vx 

    def draw(self):
        pygame.draw.rect(screen, white, (self.x, self.y, 10,5))      

##############################################################################################################################################
#declare instance of the pendulum class. Empty list for bullet objects. 
bob = pendulum(w,r,m)
bullets = []

#initial parameters before start of the game.         
b = 0.0
F = 0.0
wext = 0.05
frame = 0

#function to print text to screen
font = pygame.font.SysFont(None, 40)
def text_screen(text, colour, x, y):
    screen_text = font.render(text,True,colour)
    screen.blit(screen_text,(x,y))

def instructionscreen():
    play = True
    while play:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:  
                play = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_RETURN:
                    gameloop()
        screen.fill(black)
        text_screen('Welcome to the Pendulum Equilibrium Game!', white, 50, 75)
        text_screen('You have to maintain a score of 0!', white, 50, 125)
        text_screen('Every 5 seconds a bullet is launched.', white, 50, 175)
        text_screen('If it hits you: -20, If you dodge: +10.', white, 50, 225)
        text_screen('Change damping by Left and Right keys.', white, 50, 275)
        text_screen('Change force by Up and Down keys.', white, 50, 325)
        text_screen('Use phyiscs to escape the bullet!', white, 50, 375)
        text_screen('Feel ready? Hit Enter', white, 50, 425)
        text_screen('Get back to this screen by hitting Enter', red, 30, 475)
        text_screen('IMPORTANT: The goal is to maintain a score of 0!', red, 30, 525)

        pygame.display.update()
        clock.tick(fps)

    pygame.quit()
    quit()


def gameloop():
    play = True
    #can avoid using global variables if these are parameters in the pendulum class. but oh well....
    global b
    global F
    global frame
    while play:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:  
                play = False

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_RIGHT:
                    if b<=9.0:
                        b = b + 1.0
                if event.key == pygame.K_LEFT:
                    if b>=1.0:
                        b = b - 1.0    
                if event.key == pygame.K_UP:
                    F = F + 0.1
                    if F>=1.0:
                        F = 1.0
                if event.key == pygame.K_DOWN:
                    F = F - 0.1
                    if F<=0.0:
                        F = 0.0
                if event.key == pygame.K_RETURN:
                    instructionscreen()
                '''
                if event.key == pygame.K_x:
                    wext = wext + 0.1
                if event.key == pygame.K_z:
                    wext = wext - 0.1 
                '''

        screen.fill(black)

        #move the pendulum and the bullets. check if the bullets have left the boundary               
        bob.RK4(b,F,wext)
        
        for bull in bullets:
            bull.move()
        
        #check for pendulum-bullet collisions + destroy them + update the score
        l = 0
        while l<len(bullets):
            if (abs(bullets[l].x + 5 - screen_width/2 - bob.x)<= bob.radius + 5)  & (abs(bullets[l].y + 2.5 - screen_height/2 - bob.y)<=bob.radius + 2.5):
                bullets.pop(l)
                bob.score = bob.score - 20 
            l = l+1
        l = 0
        
        # add bullet objects to the bullet list
        if frame % 300 == 0 and frame>=300:
            bullet_x = np.random.randint(0,2)
            if bullet_x  == 0:
                bullet_x = np.random.randint(25,50)
                bullet_vx = 5.0
            elif bullet_x == 1:
                bullet_x = np.random.randint(screen_width-50, screen_width-25)
                bullet_vx = -5.0
            bullet_y = np.random.randint(screen_height/2 + bob.y-26 ,screen_height/2 + bob.y + 25)
            bullets.append(bullet(bullet_x, bullet_y, bullet_vx))

        #drawing the boundary lines
        pygame.draw.line(screen, green, (0,screen_height/4 - bob.radius), (screen_width,screen_height/4 - bob.radius))
        pygame.draw.line(screen, green, (0,3*screen_height/4 + bob.radius), (screen_width,3*screen_height/4 + bob.radius))

        
        #Drawing the pendulum
        pygame.draw.line(screen, red, (screen_width/2, screen_height/2), (screen_width/2 + bob.x,screen_height/2 + bob.y))
        bob.draw()
        pygame.draw.circle(screen, blue, (screen_width/2, screen_height/2),5)

        #Drawing the bullets
        for bull in bullets:
            bull.draw()
        
        #Destroying the bullets if they leave the screen + update the score
        l = 0
        while l<len(bullets):
            if bullets[l].x <= 0 or bullets[l].x>=screen_width:
                bullets.pop(l)
                bob.score = bob.score + 10  
            l = l+1
        l = 0

        
        #Display Text to screen
        text_screen('Time: '+ str(int(bob.t)), white, 50, 75)
        text_screen('Energy: '+ str("%.3f"%bob.E), white, 250, 75)
        text_screen('Natural w: '+ str("%.3f"%w), white, 550, 75)
        text_screen('Damping b: '+ str("%.2f"%b), white, 50, 525)
        text_screen('LEFT/RIGHT', red, 50, 565)
        text_screen('Amplitude F: '+ str("%.2f"%F), white, 300, 525)
        text_screen('UP/DOWN', red, 300, 565)
        text_screen('Score: '+ str(bob.score), white, 600, 525)
        text_screen('GOAL: 0', red, 600, 565)
        
        frame = frame + 1
        pygame.display.update()
        clock.tick(fps)
        
    pygame.quit()
    quit()

instructionscreen()
#gameloop() since my instruction screen already calls on the gameloop, I don't need to call it separately.
