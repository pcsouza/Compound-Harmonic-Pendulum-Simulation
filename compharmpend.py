#This program produces a simulation of a double pendulum with the two elements attached by springs. It uses the 4th-order dynamic step size runge-kutta approximation method with excellent results. Testing shows that even with maximum allowed error of .01 joules, the cumulative error in energy is a fraction of a percent. The user can specify custom parameters in the prompt or opt for default parameters.  

from math import pi, sin, cos, sqrt
from numpy import array, amax, zeros, sum, empty, shape
from visual import sphere, cylinder, rate, label, display
from time import sleep
import csv
import sys


# Define the primary global constants and initial parameters of the simulation using a csv datafile. Each row is represents a pendulum element with the order of parameters being:
#x, y position
#x, y velocity
#preceding spring rest length
#spring constant
#mass

def constants(filename):
	global g, m1, m2, k1, k2, l00, l01, c, targetacc, mass, spring_k, gravity, l_0
	g = 9.81   		#gravity makes the world go round
	gravity = 9.81
	c = 0 	#A constant in a small corrective force that prohibits negative
			#extension of springs
	targetacc = .001
	
	#read files
	with open(filename, 'r') as f:
		try:
			dialect = csv.Sniffer().sniff(f.read())
			f.seek(0)
			hasheader = csv.Sniffer().has_header(f.read())
			f.seek(0)
		except:
			print "Unable to interpret CSV dialect"
				
		#Count the number of pendulum elements
		for i, j in enumerate(f):
			pass
		f.seek(0)
		reader = csv.reader(f, dialect)
		if hasheader == True:
			numelements = i
			reader.next()
		else:
			numelements = i+1
		
		#Assign parameters
		r = empty((numelements, 2, 2))
		mass, l_0, spring_k = [], [], []		
		try:
			for i in range(numelements):
				line = reader.next()
				r[i, 0, 0], r[i, 0, 1] = float(line[0]), float(line[1])
				r[i, 1, 0], r[i, 1, 1] = float(line[2]), float(line[3])
				l_0.append(float(line[4]))
				spring_k.append(float(line[5]))
				mass.append(float(line[6])) 
		except:
			print "Incorrectly-formatted file"
	return r

#Function for the distance between elements. This is convenient for other calculations 
#in the simulation
def length(r):		#Calculate the distance between pendulum elements
	l = [sqrt((r[0,0,0]**2)+(r[0,0,1]**2))]
	for i in range(1, shape(r)[0]):
		l.append(sqrt(((r[i-1,0,0]-r[i,0,0])**2) + ((r[i-1,0,1]-r[i,0,1])**2)))
	return l

#dl/dx for permutations of lengths and parameters. Again, convenient for producing 
#other equations in the simulation.
def difflength(r):	#Calculate differential components of length
	dlength = empty([shape(r)[0], 2, 2])
	l = length(r)
	dlength[0, 0, :] = -r[0,0,:]/l[0]
	for i in range(1, shape(r)[0]):
		dlength[i, 0, :] = (r[i-1, 0, :] - r[i,0, :])/l[i]
	for i in range(shape(r)[0]):
		dlength[i, 1, :] = -dlength[i, 0, :]
	
	return dlength 	#indices organized by [spring element, with respect to pendulum
					#element, in the x or y coordinate]

def EqsOfMotion(r):
	l = length(r)
	dl = difflength(r)
	Eqs = empty([shape(r)[0],2,2])
	L = shape(r)[0]-1
	
	for i in range(L):
		#Motion due to springs
		Eqs[i, 1, :] = ((spring_k[i]*(l_0[i] - l[i])*(dl[i,1,:]) + spring_k[i+1]* \ 
		(l_0[i+1] - l[i+1])*(dl[i+1,0,:])) + (c*dl[i,1,:]/(l[i]**2) + \
		c*dl[i+1,0,:]/(l[i+1]**2)))/mass[i]
		#Plus gravity
		Eqs[i, 1, 1] -= gravity
	
	Eqs[L, 1, :] = ((spring_k[L]*(l_0[L] - l[L])*(dl[L,1,:])) + (c*dl[L,1,:]/(l[L]**2)))/mass[L]
	Eqs[L, 1, 1] -= gravity
	
	Eqs[:, 0, :] = r[:, 1, :]

	return Eqs

# The 4th order Runge-Kutta approximation method. Produces each individual step based on the equations of motion	
def Runge(r, step):
	k1 = step*EqsOfMotion(r)
	k2 = step*EqsOfMotion(r + .5*k1)
	k3 = step*EqsOfMotion(r + .5*k2)
	k4 = step*EqsOfMotion(r + k3)
	return (k1 + 2*k2 + 2*k3 + k4)/6

#Initializes and redraws the animation elements. The size of the spheres and rods is based off of the neutral spring lengths. Produces a text display of the percent error of net energy, which is the best metric of simulation energy because analytical equations of energy exist for the system and energy is conserved.
def Draw(r, init):
	Error = ((NetEnergy(r) - Energy0)/Energy0)
	if init == 1:	#Initializes the display window
		global pend0, pend1, pend2, rod1, rod2, ErrDisplay, scene, pend, rod
		scene = display(width = 800, height = 800)
		pend = [sphere(pos=(0,0,0), radius = .1*(l_0[0] + l_0[1]))]
		rod = [cylinder(pos = (r[0,0,0],r[0,0,1],0), axis = (-r[0,0,0], r[0,0,1], 0), radius = .025*(l_0[0] + l_0[1]))]
		for i in range(shape(r)[0]):
			pend.append(sphere(pos = (r[i,0,0],r[i,0,1],0), radius = .1*(l_0[0] + l_0[1])))
		for i in range(1, shape(r)[0]):
			rod.append(cylinder(pos = (r[i,0,0],r[i,0,1],0), axis = (r[i-1,0,0]-r[i,0,0], r[i-1,0,1]-r[i,0,1],0), radius = .025*(l_0[0] + l_0[1])))
		ErrDisplay = label(pos = (0, l_0[0] + l_0[1] + 2,0), text = "Error Percent: %f " % Error)
	else:	#Redraw for steps
		for i in range(shape(r)[0]):
			pend[i+1].pos = (r[i,0,0],r[i,0,1],0)
		rod[0].axis = (-r[0,0,0], -r[0,0,1],0); rod[0].pos = (r[0,0,0], r[0,0,1], 0)
		for i in range(1, shape(r)[0]):
			rod[i].pos = (r[i,0,0], r[i,0,1], 0)
			rod[i].axis = (r[i-1,0,0]-r[i,0,0], r[i-1,0,1]-r[i,0,1],0)
		ErrDisplay.text = "Error Percent: %f " % Error

#Calculates the maximum 5th-order error for each step, then applies an approximate error propagation formula to convert the various uncertainties into a single uncertainty for energy. The ratio of the target accuracy and the maximum error is used to adjust the step size of the simulation.
def AccuracyRatio(r, step, targetacc):
 	rtest0 = r + Runge(r, step)
	rtest1 = rtest0 + Runge(rtest0, step)
	rtest2 = r + Runge(r, 2*step)
	rerror = rtest2-rtest1
	rerror /= 30
 	dE = dEnergy(r)
 	Eprop = float(0)
 	for i in range(2):
 		for j in range(2):
 			for k in range(2):
 				Eprop += (rerror[i,j,k]*dE[i,j,k])**2
 	Eerror = sqrt(Eprop)
 	if Eerror == 0:
 		Eerror += .00001
	return targetacc/Eerror

#Net energy of any set of parameters. Convenient for error calculation and debugging.
def NetEnergy(r):
	l = length(r)
	E = 0
	for i in range(len(l)):
		E += .5*mass[i]*(r[i,1,0]**2 + r[i,1,1]**2) + .5*spring_k[i]*((l_0[i]-l[i])**2) \
		 + c/l[i] + mass[i]*gravity*r[i,0,1]

	return E

#Energy with respect to the system state. Convenient for calculating error propagation
def dEnergy(r):
	l = length(r)
	dl = difflength(r)
	dE = empty(shape(r))
	L = shape(r)[0] - 1
	for i in range(L):
		dE[i, 0, :] = -spring_k[i]*(l_0[i]-l[i])*dl[i,1,:] - spring_k[i+1]*(l_0[i] - l[i])*dl[i+1,0,:] - c*dl[i,1,:]/(l[i]**2) - c*dl[i+1,0,:]/(l[i+1]**2)
		dE[i, 0, 1] += mass[i]*gravity
	dE[L, 0, :] = -spring_k[L]*(l_0[L]-l[L])*dl[L,1,:] - c*dl[L,1,:]/(l[L]**2)
	dE[L, 0 , 1] += mass[L]*gravity
	for i in range(L):
		dE[i, 1, :] = mass[i]*r[i, 1, :]
	return dE


#Main
r = constants(sys.argv[1])		#Establish initial conditions
Energy0 = NetEnergy(r)	#Establish a global initial energy
time = 0
step = .01	#Initial step size
Draw(r, 1)	#Initialize the animation elements
runinter = 0
while runinter < 10000:
	Accuracy = AccuracyRatio(r, step, targetacc)	#Calculate the step error
	if Accuracy < 1:								#If the error is too great
 		step*=(1*Accuracy**.25)						#reduce step size and retry
 		
 	else:											#Else, use the new step and 
		r += Runge(r, step)							#Increase the step size to be
 		if step < .02:								#more efficient. Maximum stepsize
 			step *= (1*Accuracy**.25)				#imposed so the animation runs 
 		else:										#smoothly
 			step = .02
		Draw(r,0)
 		time += step
 		rate(1/step)
