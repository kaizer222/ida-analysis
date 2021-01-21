import numpy as np
class Spring:
    def __init__(self, constant, length, series_qty, parallel_qty, init_xpos, init_ypos, init_cmprs):
        self.constant=constant
        self.length=length
        self.series_qty=series_qty
        self.parallel_qty=parallel_qty
        self.init_xpos=init_xpos
        self.init_ypos=init_ypos
        self.init_cmprs=init_cmprs

    def length_total(self):
        return 0.05*self.series_qty

class Ball:
    def __init__(self, mass_ball, mass_carrier, diameter):
        self.mass_ball=mass_ball
        self.mass_carrier=mass_carrier
        self.diameter=diameter

    def totalmass(self):
        return self.mass_ball+self.mass_carrier


class LaunchConditions:
    def __init__(self, angle_in_deg, g_accel, frict_coef, init_vel):
        self.angle_in_deg=angle_in_deg
        self.g_accel=g_accel
        self.frict_coef=frict_coef
        self.init_vel=init_vel

    def angle_in_rad(self):
        return self.angle_in_deg/180*np.pi

class Air:
    def __init__(self, density, dynamic_visc):
        self.density=density
        self.dynamic_visc=dynamic_visc


Spring=Spring(1000, 0.05,2,3,0.2,0.1,4)
Ball=Ball(0.050,0.025, 0.040)
LaunchConditions=LaunchConditions(50,9.81,0.5, 0.0)
Air=Air(1.204, 1.825e-5)
