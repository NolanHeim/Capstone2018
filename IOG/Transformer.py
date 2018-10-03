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
        self.Seconds_Per_Day = 86400.0
        self.Seconds_Per_Day_Stars = 86164.0
        self.Seconds_Per_Hour = 3600.0
   
    
    def geo_2_eci(self, lat, lon, alt, JDtime): #theta in radians
        theta = self.get_theta(JDtime) 
        
        phi = theta + math.radians(lon)
        phi_dot = 2.0*math.pi/self.Seconds_Per_Day_Stars
        #phi_dot = 2.0*math.pi/self.Seconds_Per_Hour
        
        radius = alt + self.r_Earth*math.cos(math.radians(lat))
        x = radius*math.cos(phi)
        y = radius*math.sin(phi)
        z = (alt + self.r_Earth)*math.sin(math.radians(lat))
        
        vx = (-1.0)*radius*math.sin(phi)*phi_dot
        vy = radius*math.cos(phi)*phi_dot
        vz = 0.0        
        
        return [x, y, z, vx, vy, vz]
        

    def ecef_2_eci(self, x, y, z, vx, vy, vz, time, JDtime):
        theta = self.get_theta(JDtime + self.seconds_2_days(time))

        #positional arguemnts        
        xECI = x*math.cos(theta) - y*math.sin(theta)
        yECI = x*math.sin(theta) + y*math.cos(theta)
        zECI = z
        
        #veclocity arguments
        vxECI = vx*math.cos(theta) - vy*math.sin(theta)
        vyECI = vx*math.sin(theta) + vy*math.cos(theta)
        vzECI = vz
        
        return [xECI, yECI, zECI, vxECI, vyECI, vzECI]        
        

    def construct_site_matrix(self, lat, lon, alt, times, JDtime):
        matrix = []        
        
        for i in range(0,len(times)):
            JD_new = JDtime + self.seconds_2_days(times[i])
            
            siteInfoLine = self.geo_2_eci(lat, lon, alt, JD_new)            
            
            matrix.append(siteInfoLine)
        
        return matrix
    
    
    def seconds_2_days(self, t_sec):
        return t_sec / self.Seconds_Per_Day
    
    
    def get_theta(self, JDtime):
        day = JDtime - 2451545.0
        GMST = 18.697374558 + 24.06570982441908*day
#        return GMST*(2.0*math.pi/24.0)
        return GMST*(2.0*math.pi/self.Seconds_Per_Day)
