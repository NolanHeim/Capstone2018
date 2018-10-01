# -*- coding: utf-8 -*-
#
# Transformer.py
#
# Initial Creation Date: 09/30/2018
#
# Written by Jordan Jones and Nolan Heim
#

import math

class Transformer:
    
    def __init__(self):
        self.r_Earth = 6378137.0 #m
   
    
    def geo_2_eci(self, lat, lon, alt, JDtime): #theta in radians
        D = JDtime - 2451545.0
        GMST = 18.697374558 + 24.06570982441908*D
        theta = GMST*(2*math.pi/24.0) 
        
        phi = theta + math.radians(lon)
        
        radius = alt + self.r_Earth*math.cos(math.radians(lat))
        x = radius*math.cos(phi)
        y = radius*math.sin(phi)
        z = (alt + self.r_Earth)*math.sin(math.radians(lat))
        
        return [x, y, z]
        
    def ecef_2_eci(self, x, y, z, vx, vy, vz, JDtime):
        
        D = JDtime - 2451545.0
        GMST = 18.697374558 + 24.06570982441908*D
        theta = GMST*(2*math.pi/24.0)

        #positional arguemnts        
        xECI = x*math.cos(theta) - y*math.sin(theta)
        yECI = x*math.sin(theta) + y*math.cos(theta)
        zECI = z
        
        #veclocity arguments
        vxECI = vx*math.cos(theta) - vy*math.sin(theta)
        vyECI = vx*math.sin(theta) + vy*math.cos(theta)
        vzECI = vz
        
        return [xECI, yECI, zECI, vxECI, vyECI, vzECI]        
        

