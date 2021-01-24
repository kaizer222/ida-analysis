#!/usr/bin/env python
# coding: utf-8

# In[2]:


# import all modules required
import numpy as np              # for the arrays  
import matplotlib.pyplot as plt # for the plotting
import scipy.integrate as si    # to integrate numerically the equations of motion
import scipy.optimize as so     # for the bisection algorithm
import pandas as pd

# change to qt for plots in a new window
#get_ipython().run_line_magic('matplotlib', 'qt')


# ### Parameters (spring)

# In[ ]:


# mass of squash ball
m = 0.025 # kg

# diameter of squash ball
d =  0.040 # mm

# mass of carrier
mc = 0.025 # kg

# mass of ball + carrier
mt = m + mc

# Number of springs in series
Ns = 1

# Number of springs in parallel
Np = 1

# Length of spring
L0 = 0.1016 * Ns # m
# Spring constant
k = 794/Ns*Np # N/m
# maximumum compression distance = 0.082 m

# Launch angle
theta = 67/180*np.pi # rad

# Gravitational acceleration
g = 9.81 # m/s^2

# Friction factor
f = 0.5 # N/(m/s)

# Spring position
xs = 0.2 # m
ys = 0.1 # m


# ### Parameters (Trajectory)

# In[ ]:


# density and dynamic viscosity of air
rho = 1.204    # kg/m^3
mu  = 1.825e-5 # Pa s 


# ### Initial Conditions for spring

# In[ ]:


# Initial spring compression
s0 = 0.02 # m

# Initial velocity
v0 = 0.0 # m/s

y0 = [s0, v0]


# ### Spring equation

# In[ ]:


# Spring equation
def spring(y, t):
    # unpack
    s, v = y
    # dsdt
    dsdt = -v
    # dvdt
    dvdt = k*s/mt - g*np.sin(theta) - f*v/mt
    
    return dsdt, dvdt


# ### Find maximum launching velocity

# In[ ]:


def maxVelocity(t_end):
    t = np.linspace(0, t_end, 100)
    sol = si.odeint(spring, y0, t)
    v, a = spring(sol[-1], t_end)
    return a


# In[ ]:


plt.clf()
y = []
for s0 in np.array([0.045, 0.060, 0.075])*Ns:
    y0 = [s0, v0]

    t_end = so.newton(maxVelocity, 0.013)

    t = np.linspace(0,t_end,100)

    sol = si.odeint(spring, y0, t)

    s = sol[:,0]
    v = sol[:,1]
    y.append(v[-1])
    plt.plot(t, v, label='s0=%0.2f m' % s0)
    
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('Velocity (m/s)')

#print(y)


# ### Define functions/equations that needed for trajectory

# In[ ]:


# obtain drag coefficient from Reynolds number
def get_CD(Re):
    return 24/Re * (1 + 0.27*Re)**0.43 + 0.47 * (1 - np.exp(-0.04*Re**0.38))


# In[ ]:


# compute right hand side of governing equations at a generic
# position u, defining the x, y, Vx and Vy state variables
def squashball(u, t):
    # unpack u
    x, y, Vx, Vy = u
    
    # find theta
    theta = np.arctan2(Vy, Vx)
    
    # magnitude of velocity
    V_mag = np.sqrt(Vx**2 + Vy**2)
    
    # compute Reynolds number
    Re = rho * V_mag * d / mu
    
    # drag coefficient
    CD = get_CD(Re)
    
    # calculate drag (AND HERE IS WHERE THE MISTAKE WAS!!!)
    D_mag = 0.5 * rho * V_mag**2 * CD * d**2
    
    # calculate drag components
    Dx = D_mag * np.cos(theta)
    Dy = D_mag * np.sin(theta)
    
    # calculate acceleration
    ax = - Dx/m
    ay = - Dy/m - g
    
    # return velocity and acceleration
    return [Vx, Vy, ax, ay]


# ### Trajectory 

# In[ ]:


# return t, x, y, Vx, Vy given final integral time T
def trajectory(x0, y0, Vx0, Vy0, T=1):
    # pack IC into u0
    u0 = [x0, y0, Vx0, Vy0]
    
    # define array of times
    t  = np.linspace(0, T, 101)
    
    # integrate forward
    out = si.odeint(squashball, u0, t)

    # unpack data
    x  = out[:, 0]
    y  = out[:, 1]
    Vx = out[:, 2]
    Vy = out[:, 3]
    
    return t, x, y, Vx, Vy


# ### Trajectory with stop 

# In[ ]:


# return t, x, y, Vx, Vy from given IC with y[-1] = 0
def trajectory_with_stop(x0, y0, Vx0, Vy0):

    # return y at t=T from give initial conditions
    def find_y_at_T(T):
        t, x, y, Vx, Vy = trajectory(x0, y0, Vx0, Vy0, T=T)
        return y[-1]
    
    # find T_stop such that y(T_stop) = 0
    # find T_stop such that find_y_at_T(T_stop) = 0
    T_stop = so.newton(find_y_at_T, 2)
    
    return trajectory(x0, y0, Vx0, Vy0, T=T_stop)


# ### Define functions for calculating distance

# In[ ]:


# compute velocity components from velocity magnitude and angle
def mag2comp(V_mag, theta):
    return V_mag * np.cos(theta), V_mag * np.sin(theta)


# In[ ]:


# return distance given initial velocity (at 45 degrees)
def distance(V_mag, theta):
    # find velocity components
    Vx, Vy = mag2comp(V_mag, np.deg2rad(theta))

    # solve forward
    t, x, y, Vx, Vy = trajectory_with_stop(0, 0, Vx, Vy)
    
    # return final distance
    return x[-1]


# In[3]:


# initial velocities
V_mags = y
angle = np.linspace(0, 90, 91)

xtotal=[]
for v in V_mags:
    print ('Started outer loop',v)
    x=[]
    for theta in angle:
        c = distance(v, theta)
        x.append(c)
    xtotal.append(x)        
'''
x1 = x[0:10]
x2 = x[10:20]
x3 = x[20:30]
x4 = x[30:40]
x5 = x[40:50]
x6 = x[50:60]
x7 = x[60:70]
x8 = x[70:80]
x9 = x[80:90]
x10 = x[90:100]

xtotal = []
xtotal.append(x1)
xtotal.append(x2)
xtotal.append(x3)
xtotal.append(x4)
xtotal.append(x5)
xtotal.append(x6)
xtotal.append(x7)
xtotal.append(x8)
xtotal.append(x9)
xtotal.append(x10)
'''
    #converts x4,x7,x10 lists as Pandas dataframe, and save as csv.

output_dict={'Angle':angle,'45mm':xtotal[0],'60mm':xtotal[1],'75mm':xtotal[2]}
output_dict=pd.DataFrame(output_dict)
output_dict.to_csv('outputtable.csv')
# In[4]:


plt.clf()
plt.figure(4)

plt.plot(angle, xtotal[0], label='v(45mm) = %0.2f m/s' % V_mags[0])

plt.plot(angle, xtotal[1], label='v(60mm) = %0.2f m/s' % V_mags[1])

plt.plot(angle, xtotal[2], label='v(75mm) = %0.2f m/s' % V_mags[2])


plt.legend(loc='best')
plt.grid(1)
plt.xlabel("angle [degrees]")
plt.ylabel("distance [m]")
plt.show()

# In[ ]:




