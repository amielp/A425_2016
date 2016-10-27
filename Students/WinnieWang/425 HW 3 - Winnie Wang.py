import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import scipy.integrate as integrate

print "425 HW 3 - Hubble's Constant"
print "Winnie Wang, Prof. Andrew Connelly"
print
print "This is Question 2."
#Constants:
H_0 = 2.268e-18 #Hubble constant scaled to 1/s
H_02 = 0.07154 #Hubble constant scaled to 1/Gyr
c = 2.9979e8 #Speed of light in m/s
scale = 3.241e-23 #Scaling factor of meters to Mpc

size = 1000 #arbitrary number of steps to evaluate z
z = np.linspace(0,10,size) #linear space of z to integrate against

#H(z) solved:
def Hubble(z,r,m,l,o,H):
    return (H)*((r*((z+1)**4)+m*((z+1)**3)+(1-o)*((z+1)**2)+l))**(0.5)
    
#Part A: Where 0_r = 0, 0_m = 0.3, 0_l = 0.7
#0_o = 0 for a flat universe:
HubblePartA = lambda z: c/(Hubble(z,0,0.3,0.7,0,H_0)) #For Comoving Distance
resultPartA = np.array([])

Hubble2PartA = lambda z: 1/(Hubble(z,0,0.3,0.7,0,H_02)*(1+z)) #For Age
result2PartA = np.array([])

#Part B: Where 0_r = 0, 0_m = 10, 0_l = 0 
#0_o = 1 for a closed universe:
HubblePartB = lambda z: c/(Hubble(z,0,10,0,1,H_0)) #For Comoving Distance
resultPartB = np.array([])

Hubble2PartB = lambda z: 1/(Hubble(z,0,10,0,1,H_02)*(1+z)) #For Age
result2PartB = np.array([])

for i in z:
    integrate_A = integrate.quad(HubblePartA, 0, i)
    integrate_A2 = integrate.quad(Hubble2PartA, 0, i)
    
    resultPartA = np.append(resultPartA, integrate_A[0]*scale)
    result2PartA = np.append(result2PartA, integrate_A2[0])
    
    integrate_B = integrate.quad(HubblePartB, 0, i)
    integrate_B2 = integrate.quad(Hubble2PartB, 0, i)
    
    resultPartB = np.append(resultPartB, integrate_B[0]*scale)
    result2PartB = np.append(result2PartB, integrate_B2[0])
    
#Check to see if it works:    
#print resultPartA
#print resultPartB
#print result2PartA
#print result2PartB
    
print "Part A:"
plt.plot(z, resultPartA, "g--")
plt.plot(z, resultPartB, "b--")
plt.title("Redshift vs. Comoving Distance")
plt.xlabel("Redshift, z")
plt.ylabel("Comoving Distance in Mpc")
plt.show()

print "Part B:"
plt.plot(z, result2PartA, "g-")
plt.plot(z, result2PartB, "b-")
plt.title("Redshift vs. Age")
plt.xlabel("Redshift, z")
plt.ylabel("Age in Giga-years")
plt.show()

print "This is Question 3."

#This is the evaluated age of (benchmark) flat universe from Part A:
#With 0_m = 0.3, and 0_l = 0.7
benchmark_Age_A = integrate.quad(Hubble2PartA, 0, size)[0]

#To solve for 0_r, set 0_r = 0 in the integral for age, then evaluate with a conditional loop:
r2 = 0.0 #radiation density; given as zero

Hubble3 = lambda z: 1.0/(Hubble(z,r2,0,0,-1.0,H_02)*(1.0+z)) #For Age
Hubble3_Age = integrate.quad(Hubble3, 0, size)[0]

while(Hubble3_Age >= benchmark_Age_A):
    r2 += 0.0001 #step size; has a resolution of 1e-4 because larger resolutions won't work
    Hubble3 = Hubble3 = lambda z: 1.0/(Hubble(z,r2,0,0,-1.0,H_02)*(1.0+z))
    Hubble3_Age = integrate.quad(Hubble3, 0, size)[0]
    
print "This is the benchmark age obtained with O_m = 0.3 and O_lambda = 0.7: ", benchmark_Age_A, "Gyr."
print "This is the age of matter only open universe: ", Hubble3_Age, "Gyr."
print "Solved O_r for the open universe: ", r2
    
    
"""   
Attempted to integrate through the Newton Method:
def deriHubble(x,r,m,l,naught):
    return 70*((4*r*(x**3)+3*m*(x**2)+2*(1-naught)*(x)))/((r*(x**4)+m*(x**3)+(1-naught)*(x**2)+l))**(0.5)
    
def HubbleSolver(func,ini,deri,it):
    return optimize.newton(func,ini,fprime=deri,maxiter=it) 

#Part A:    
Hubble1 = Hubble(x=x,r=0,m=0.3,l=0.7,naught=1)
deriHubble1 = deriHubble(x=x,r=0,m=0.3,l=0.7,naught=1)
    
HubbleNewtonFunc = HubbleSolver(Hubble1,1,deriHubble1,10)
print "Part A: Root obtained from packaged Newton Method: ", HubbleNewtonFunc
"""
    
