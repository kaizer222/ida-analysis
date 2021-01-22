import param as prm
import numpy as np
import matplotlib.pyplot as plt

# V --> Re --> Cd --> D
# Vinit, Ainit, density and viscousity given
# need: xpos, ypos rel to gorund
# we want total distance
# take into account height difference

#total time to reach a specific height should give us the horizontal distance.
#max height when Vy=0, when angle is zero

muzzle_vel=[16,6,4,2] #for ezample
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

        
        
