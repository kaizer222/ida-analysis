import param as prm
import numpy as np


def spring_eqn(cmprs_len):#take in compressed length, output muzzle vel. and accel

    

    muzzle_vel=-v

    spring_force=prm.Spring.constant*cmprs_len
    g_force=prm.LaunchConditions.g_accel*np.sin(prm.LaunchConditions.angle_in_rad())
    frct_force=prm.LaunchConditions.frict_coef*v

    muzzle_accel=(spring_force-g_force-frict_force)/prm.Ball.totalmass()

    return muzzle_vel, muzzle_accel

def tube_vel_eqn(): #equation for velocity in launch tube

    
    return 

def plot_func():
    
