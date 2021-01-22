import param as prm
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

    -Ball Cd:
        use prm.Ball.cd

    -Ball reference area:
        use prm.Ball.ref_area()

    -Ball initial x coord:
        use prm.Ball.init_xpos

    -Ball initial y coord:
        use prm.Ball.init_ypos

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

    -Height of target relative to ball initial position
        use prm.LaunchConditions.target_ypos

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
muzzle_vel=[]
plt.clf()
for s0 in np.array([0.05, 0.06, 0.07, 0.08])*prm.Spring.series_qty:
    y0 = [s0, prm.LaunchConditions.init_vel]

    t_end = so.newton(maxVelocity, 0.013)

    t = np.linspace(0,t_end,100)

    sol = si.odeint(spring_eqn, y0, t)

    s = sol[:,0]
    v = sol[:,1]
    muzzle_vel.append(v[-1])
    plt.plot(t, v, label='s0=%0.2f m' % s0)

print (muzzle_vel)
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('Velocity (m/s)')
plt.show()

#---------------------------------Trajec---------------------------------#
#muzzle_vel=[16,6,4,2] #for ezample
angle=np.linspace(10,80,15)*(np.pi/180)

def drag(v, angle):
    D=0.5*prm.Air.density*(v**2)*prm.Ball.ref_area()*prm.Ball.cd
    Dx=D*np.cos(angle)
    Dy=D*np.sin(angle)
    
    return D, Dx, Dy

def time_up(V_init, angle):
    D,Dx,Dy = drag(V_init, angle)
    t1=-1*(V_init*np.sin(angle))/(-1*prm.LaunchConditions.g_accel-Dy/prm.Ball.totalmass())    
    return t1


def maxheight(V_init,angle):
    D, Dx, Dy=drag(V_init,angle)
    s=-1*(V_init*np.sin(angle)**2)/(2*prm.LaunchConditions.g_accel-Dy/prm.Ball.totalmass())
    
    return s

def time_down(V_init,angle):
    #end height+max height
    D,Dx,Dy=drag(V_init,angle)
    total_height=(prm.Ball.init_ypos-
                  prm.LaunchConditions.target_ypos)+maxheight(V_init,angle)
    t2=(2*total_height)/(prm.LaunchConditions.g_accel-Dy/prm.Ball.totalmass())**0.5
    return t2

def max_distance(V_init,angle):
    D,Dx,Dy=drag(V_init,angle)
    t=time_up(V_init, angle)+time_down(V_init, angle)
    s=(V_init*np.cos(angle))*t+0.5*(Dx/prm.Ball.totalmass())*(t**2)
    return s

#distance against angle
def show_dva():
    for i in muzzle_vel:
        d=[]
        for a in angle:
            d.append(max_distance(i,a))
        plt.plot(angle,d)
        print(d)
    plt.xlabel('Angle(rad)')
    plt.ylabel('Max Distance')
    plt.show()

show_dva()
