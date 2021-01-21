import param as prm
import spring, trajectory
import matplotlib.pyplot as plt
import scipy.integrate as si
import scipy.optimize as so
import numpy as np


#--------------------------------README----------------------------------#
'''
List of parameters:

Spring:
    -Spring constant
        use prm.Spring.constant
        
    -One spring length only
        use prm.Spring.length
        
    -Spring length considering series spring (function)
        use prm.Spring.length_total()

    -Number of series spring
        use prm.Spring.series_qty

    -Number of parallel spring
        use prm.Spring.parallel_qty

    -Initial x position of spring
        use prm.Spring.init_xpos

    -Initial y position of spring
        use prm.Spring.init_ypos

    -Initial spring compression:
        use prm.Spring.init_cmprs

Ball:
    -Ball mass
        use prm.Ball.mass_ball

    -Carrier mass
        use prm.Ball.mass_carrier

    -Ball diameter
        use prm.Ball.diameter

    -Total mass (function)
        use prm.Ball.totalmass()

Launch Conditions (refered to as Initial Conditions in Template):
    -Launch angle (in degrees):
        use prm.LaunchConditions.angle_in_deg

    -Launch angle (in radians) (function):
        use prm.LaunchConditions.angle_in_rad()

    -Gravitational accleration
        use prm.LaunchConditions.g_accel

    -Friction factor
        use prm.LaunchConditions.frict_coef

    -Initial velocity
        use prm.LaunchConditions.init_vel

Properties of air:
    -Density:
        use prm.Air.density

    -Dynamic Viscousity
        use prm.Air.dynamic_visc

ALL PARAMETERS EDITABLE IN PARAM.PY
UPDATE HERE FOR ADDITIONAL PARAMETERS
'''
#------------------------------------------------------------------------#

#---------------------------------SPRING---------------------------------#

#SPRING EQUATION
y0=[prm.Spring.init_cmprs, prm.LaunchConditions.init_vel]
def spring_eqn(y, t):
    s, v = y
    # dsdt
    dsdt = -v
    
    # dvdt
    k=prm.Spring.constant
    mt=prm.Ball.totalmass()
    g=prm.LaunchConditions.g_accel
    theta=prm.LaunchConditions.angle_in_rad()
    f=prm.LaunchConditions.frict_coef
    v=prm.LaunchConditions.init_vel
    
    dvdt = k*s/mt - g*np.sin(theta) - f*v/mt
    
    return dsdt, dvdt

def maxVelocity(t_end):
    t = np.linspace(0, t_end, 100)
    sol = si.odeint(spring_eqn, y0, t)
    v, a = spring_eqn(sol[-1], t_end)
    return a
     
plt.clf()
muzzle_velocity=[]
plt.clf()
for s0 in np.array([0.01, 0.02, 0.03, 0.04])*prm.Spring.series_qty:
    y0 = [s0, prm.LaunchConditions.init_vel]

    t_end = so.newton(maxVelocity, 0.013)

    t = np.linspace(0,t_end,100)

    sol = si.odeint(spring_eqn, y0, t)

    s = sol[:,0]
    v = sol[:,1]
    muzzle_velocity.append(v[-1])
    #plt.plot(t, v, label='s0=%0.2f m' % s0)
    
#plt.legend()
#plt.xlabel('time (s)')
#plt.ylabel('Velocity (m/s)')
#plt.show()

print (muzzle_velocity)

#------------------------------------------------------------------------#
def Vx(V, angle):
    return V*np.cos(angle)
def Vyi(V, angle):
    return V*np.sin(angle)

def max_distance(V, angle):
    Tmax=2*Vyi(V, angle)/prm.LaunchConditions.g_accel
    Dmax=Vx(V,angle)*Tmax
    return Dmax

for i in muzzle_velocity:
    d=[]
    angle=np.linspace(10, 80, num=15)
    for a in angle:
        d.append(max_distance(i, a))
    #print (len(d))
    #print (len(angle))
    plt.plot(angle, d)
#print (angle)
plt.show()
    

