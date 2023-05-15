# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:05:58 2020

@author: lenovo
"""

import matplotlib.pyplot as plt


#function definitions

def f1(x,y): #x coordinate for rocket's acceleration
    
    a = -(G*MassE*x)
    b = ((x**2) + (y**2))**(3/2)
    
    return a/b
    
def f2(x,y): #y coordinate for rocket's acceleration
    
    c = (G*MassE*y)
    d = ((x**2) + (y**2))**(3/2)
    
    return -c/d

def f3(x,y): #x coordinate for rocket's acceleration

    a = f1(x,y)
    b = (G*MassM*(x - distance))
    c = ((x - distance)**2 + (y - distance)**2)**(3/2)
    
    return a - b/c

def f4(x,y): #y coordinate for rocket's acceleration
    
    d = f2(x,y)
    e = (G*MassM*(y - distance))
    f = ((x - distance)**2+(y - distance)**2)**(3/2)
 
    return d - e/f

def RungeKuttaA(x,vx,y,vy,t): #next value of x, y, vx, vy
    
    xpos, ypos, TotalE, KinE, PotE, time = [],[],[], [], [], []
    
    while (t<t0):

        x1 = vx
        y1 = vy
        vx1 = f1(x,y)
        vy1 = f2(x,y)
        x2 = vx + (a*(vx1))/2
        y2 = vy + (a*(vy1))/2
        vx2 = f1(x+(a*x1)/2,y+(a*y1)/2)
        vy2 = f2(x+(a*x1)/2,y+(a*y1)/2)
        x3 = vx + (a*vx2)/2
        y3 = vy + (a*vy2)/2
        vx3 = f1(x+(a*x2)/2,y+(a*y2)/2)
        vy3 = f2(x+(a*x2)/2,y+(a*y2)/2)
        x4 = vx + (a*vx3)
        y4 = vy + (a*vy3)
        vx4 = f1(x+(a*x3),y+(a*y3))
        vy4 = f2(x+(a*x3),y+(a*y3))
        
        x += (a/6)*(x1+2*x2+2*x3+x4)
        y += (a/6)*(y1+2*y2+2*y3+y4)
        vx += (a/6)*(vx1+2*vx2+2*vx3+vx4)
        vy += (a/6)*(vy1+2*vy2+2*vy3+vy4)
        t += a
        
        KE = (1/2)*m*((vx**2)+(vy**2))
        PE = -(G*MassE*m)/(((x**2)+(y**2))**(1/2))
        E = KE + PE
        
        xpos.append(x)
        ypos.append(y)
        TotalE.append(E)
        KinE.append(KE)
        PotE.append(PE)
        time.append(t)
        
    return xpos, ypos, TotalE, KinE, PotE, time

def RungeKuttaB(x,vx,y,vy,t): #next value of x, y, vx, vy

    xpos1, ypos1, TotalE1, time1 = [],[],[], []
    
    while (t<t0):
            x1 = vx
            y1 = vy
            vx1 = f3(x,y)
            vy1 = f4(x,y)
            x2 = vx + (a*(vx1))/2
            y2 = vy + (a*(vy1))/2
            vx2 = f3(x+(a*x1)/2,y+(a*y1)/2)
            vy2 = f4(x+(a*x1)/2,y+(a*y1)/2)
            x3 = vx + (a*vx2)/2
            y3 = vy + (a*vy2)/2
            vx3 = f3(x+(a*x2)/2,y+(a*y2)/2)
            vy3 = f4(x+(a*x2)/2,y+(a*y2)/2)
            x4 = vx + (a*vx3)
            y4 = vy + (a*vy3)
            vx4 = f3(x+(a*x3),y+(a*y3))
            vy4 = f4(x+(a*x3),y+(a*y3))
            
            x += (a/6)*(x1+2*x2+2*x3+x4)
            y += (a/6)*(y1+2*y2+2*y3+y4)
            vx += (a/6)*(vx1+2*vx2+2*vx3+vx4)
            vy += (a/6)*(vy1+2*vy2+2*vy3+vy4)
            t += a
            
            KE = (1/2)*m*((vx**2)+(vy**2))
            PE = -(G*MassM*m)/((x**2)+(y**2))
            E = KE + PE
        
            xpos1.append(x)
            ypos1.append(y)
            TotalE1.append(E)
            time1.append(t)
        
    return xpos1, ypos1, TotalE1, time1

#value definitions
    
RadiusM = 1737e3 #radius of the moon
RadiusE = 6371e3 #radius of the Earth
distance = 3844e5
m = 1e4 #mass of rocket
G = 6.67e-11 #gravitational constant
MassE = 5.97e24 #mass of the Earth
MassM = 7.35e22 #mass of the Moon
vx = 11000 #initial velocity in x direction 
vy = 0 #initial velocity in y direction 
a = 10 #time step
x = 0 #initial position in the x coordinate
y = -6371e3 #initial position in the y coordinate
t = 0 #start at time is zero
t0 = 1E4

#menu

MyInput = '0'

while MyInput != 'q':

    MyInput = input('Enter a choice, "a", "b", or "q" to quit: ')
    print('You entered the choice ', MyInput)
    
#section A
    
    if MyInput == 'a':
        
        #define lists used in plotting
        xpos, ypos, TotalE, KinE, PotE, time = RungeKuttaA(x,vx,y,vy,t)
        
        #plot displacement y against x
        plt.plot(xpos,ypos)
        plt.title('Eliptical Orbit of the Rocket')
        plt.xlabel('Horizontal Displacement x (m)')
        plt.ylabel('Vertical Displacement y (m)')
        plt.show()
        
        #plot energy against time
        plt.plot(time,TotalE)
        plt.title('Total Energy as a Function of Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Total Energy (J)')
        plt.show()
        
        #plot kinetic energy against time
        plt.plot(time,KinE)
        plt.title('Kinetic Energy as a Function of Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Kinetic Energy (J)')
        plt.show()
        
        #plot gravitational potential energy against time
        plt.plot(time,PotE)
        plt.title('Potential Energy as a Function of Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Potential Energy (J)')
        plt.show()
        
#section B
        
    elif MyInput == 'b': 
        
        #define lists used in plotting
        xpos1, ypos1, energy1, time1 = RungeKuttaB(x,vx,y,vy,t)
        
        #plot energy against time
        plt.plot(time1,energy1)
        plt.title('Total Energy as a Function of Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Total Energy (J)')
        plt.show()
        
        #plot displacement y displacement against x displacement
        EF= plt.Circle((0,0),RadiusE)
        MF= plt.Circle((0,distance),RadiusM, color = 'r')
        orbit=plt.figure()
        ax=orbit.add_subplot(1,1,1)
        ax.add_patch(EF)
        ax.add_patch(MF)
        plt.plot(xpos1,ypos1)
        plt.title('Rocket Launch from Earth')
        plt.xlabel('Horizontal Displacement (m)')
        plt.ylabel('Vertical Displacement (m)')
        plt.plot(xpos1,ypos1)
        plt.axis('equal')
        
        plt.ylim([-6E7,4e8])
        plt.show()
        
#quitting program

    elif MyInput == 'q':
        
        print('You have chosen to quit.')
        
print('Thank you for using this program. Good bye.')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    