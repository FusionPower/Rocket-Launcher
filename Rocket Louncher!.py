#!/usr/bin/env python3
#
# n-body.py Solve the n-body problem using Newton
# 
# Copyright (C) 2019  Rocket-Launcher
#                      
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


#########################Available on Github##############################
#### This code was made with help of Victor de la Luz's code for N-Body###
##########################################################################
from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
import math

G=6.674e-11         #m^3kg^-1s^-2


class Particle:
	#Each particle will have a position, velocity and a mass. 
    def __init__(self, p, v, m, dt=1):
        self.p = p #position
        self.v = v #velocity
        self.m = m #mass
        self.dt = dt
		#trajectory and time are lists so they can be graphed easily
        self.trajectory = [p]
        self.time = [0.0]

	#to change dt
    def setdt(self,dt):
        self.dt = dt

	#Finds the radius (distance) between two particles:
	#The one called by the class and another particle p1.
    def computeR(self,p1):
        r = math.sqrt( (p1[0]-self.p[0])**2 + (p1[1]-self.p[1])**2 + (p1[2]-self.p[2])**2)
        return r

	#Calculates the unit vector that will help determine the force that
	#p1 exerts on this particle.
    def computeU(self,p1):
        u=[0,0,0]
        i=0
        for a,b in zip(self.p,p1):
            u[i] = b - a
            i+=1
        return u
	
    def getPosition(self):
        return self.p

    def getVelocity(self):
        return self.v
	
	
    def getKineticEnergy(self):
        k= (1/2)*self.m*(math.sqrt( self.v[0]^2 +self.v[1]^2+self.v[2]^2))
        return k    
	

	#Find velocity
    def computeV(self,B):
        r = self.computeR(B.p)
        u = self.computeU(B.p)


        Vx=(G*B.m*self.dt/(r**3))*u[0]
        Vy=(G*B.m*self.dt/(r**3))*u[1]
        Vz=(G*B.m*self.dt/(r**3))*u[2]

        return [Vx,Vy,Vz]


	#Helps update the velocity
    def updateV(self,w):
        self.v[0] += w[0]
        self.v[1] += w[1]
        self.v[2] += w[2]
        
	#Helps update the position using position+velocity*dt 
    def updatePosition(self,time,save):        
        self.p = [self.p[0]+ (self.v[0]) *dt,self.p[1]+ (self.v[1])*dt,self.p[2]+ (self.v[2])*dt]
        if save:
            self.time.append(time)
            self.trajectory.append(self.p)

    def getTrajectory(self):
        return self.time, self.trajectory
   
     
class Potential:
    
    def __init__(self, system, dt):
        self.system = system #set of Particles
        self.dt = dt #how often we measure changes

	#updates velocities and position in all possible pairs pof particles. Only two in this case.
    def integrate(self,time,save):
        for particle in self.system:
            for other in self.system:
                if other != particle:
                    velocity = particle.computeV(other)
                    particle.updateV(velocity)
        for particle in self.system:
            particle.updatePosition(time,save)

        return self.system



#Measures the trajectory during 10 days with a measurement every second
lenTime=3600.0*24*10  #sec
dt=1      #sec  
  

#given a degreee of launch and a starting velocity, ignoring air resistance,
#this function yields the trajectory of the rocket (and the Earth).
def test_rocket(degrees,velocity):
	print("Degrees:",degrees,"Starting velocity:", velocity,"m/s")	
	

	#obtains X component and Y component of the launch (z is ignored because Earth is a stationary sphere in this case)
	compx=velocity*math.cos(degrees*math.pi/180)
	compy=velocity*math.sin(degrees*math.pi/180)
	print("X component:",compx,"m/s", "Y component:",compy, "m/s")
	
	#particle=(position in meters, velocity in meters/s, mass in kg and dt in seconds)
	earth = Particle([0,0,0], [0, 0, 0], 6e24,dt)   
	satelite = Particle([0,6371002,0],[compx,compy,0],75e4,dt)
	
	#declare all particles in the system (open posibility of multiple systems)
	particles = [earth,satelite]


	n_steps = int(lenTime/dt)
	#declare the system
	twoBody = Potential(particles,dt)

	skip=0
	save=False
	for t in range(1,n_steps):
		#save every 1000 iterations so information doesn't overload the system
		#for more reliable set condition to skip==0 and skip+=0
		if skip == 1000: 
		    skip=0
		    save=True
		system = twoBody.integrate(float(t)*dt,save)
		save=False
		skip += 1
	
	#defining the plot figure
	fig = plt.figure()

	ax = fig.add_subplot(111, projection='3d')
	ax.set_xlim3d(-5e7, 5e7)
	ax.set_ylim3d(-5e7,5e7)
	ax.set_zlim3d(-5e7,5e7)
	
	#color of trajectory
	c=['g']
	
	#Find if the trajectory was unsuccesful 
	#i is useful for multiple colored trajectories
	i=0
	time, trajectory = particles[1].getTrajectory()
	for x,y in enumerate(trajectory):
		if np.linalg.norm(y)<6371000:
			print("The rocket crashed!")
			plt.close()
			return 1
		ax.scatter(y[0], y[1], y[2], marker='o',c=c[i])
	i=i+1



	#add a sphere for earth#
	u = np.linspace(0, 2 * np.pi, 13)
	v = np.linspace(0, np.pi, 7)

	x = 6371000 * np.outer(np.cos(u), np.sin(v))
	y = 6371000 * np.outer(np.sin(u), np.sin(v))
	z = 6371000 * np.outer(np.ones(np.size(u)), np.cos(v))

	xdata = scipy.ndimage.zoom(x, 3)
	ydata = scipy.ndimage.zoom(y, 3)
	zdata = scipy.ndimage.zoom(z, 3)

	ax.plot_surface(xdata, ydata, zdata, rstride=2, cstride=2, color='b', shade=0)

	print("Mision Accomplished!")
	plt.show()

	return 0



#This tool helps test multiple missions varying velocity and angle of launch
for velocity in range(8000,10000,500): #Specify initial velocity, maximum velocity and skips in every simulation
	for degrees in range(0,90,10):  # Specify starting angle, ending angle and skips
		x=test_rocket(degrees,velocity)
		#if x==1:    # Uncomment if you want to simulate until a mission fails
		#	break
		#if x==0:    # Uncomment if you want to simulate until a mission is successful
		#	break


