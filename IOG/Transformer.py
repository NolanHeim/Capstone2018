# -*- coding: utf-8 -*-
#
# Transformer.py
#
# Initial Creation Date: 09/30/2018
#
# Written by Jordan Jones and Nolan Heim
#

import math
import numpy as np

class Transformer:
    
    def __init__(self):
        self.r_Earth = 6378137.0 #m
        self.Seconds_Per_Day = 86400.0
        self.Seconds_Per_Day_Stars = 86164.0
        self.Seconds_Per_Hour = 3600.0
   
    
    def geo_2_eci(self, lat, lon, alt, times, JDtime): #theta in radians
        print('JD time: ' + str(JDtime))
        print(times)
        print(self.seconds_2_days(times))
        
        theta = self.get_theta(JDtime + self.seconds_2_days(times)) 
        print(theta)
        print("Length of time in Transformer is "+str(len(times)))        
        print("Length of theta is then "+str(len(theta)))
        print("Theta difference is "+str(theta[800] - theta[0]))
        
        phi = theta + np.radians(lon)
        phi_dot = 2.0*np.pi/self.Seconds_Per_Day_Stars
        
        radius = alt + self.r_Earth*math.cos(math.radians(lat))
        x = radius*np.cos(phi)
        y = radius*np.sin(phi)
        z = (alt + self.r_Earth)*np.sin(np.radians(lat))*np.ones(theta.shape)
        
        print("length of x is "+str(len(x)))        
        
        vx = (-1.0)*radius*np.sin(phi)*phi_dot
        vy = radius*np.cos(phi)*phi_dot
        vz = np.zeros(theta.shape) #Needs to be same dimension but zero

        return np.column_stack([x, y, z, vx, vy, vz])
        

    def ecef_2_eci(self, x, y, z, vx, vy, vz, time, JDtime):
        print('JD time: ' + str(JDtime))
        theta = self.get_theta(JDtime + self.seconds_2_days(time))
        
        #positional arguemnts        
        xECI = x*np.cos(theta) - y*np.sin(theta)
        yECI = x*np.sin(theta) + y*np.cos(theta)
        zECI = z
        
        #veclocity arguments
        vxECI = vx*np.cos(theta) - vy*np.sin(theta)
        vyECI = vx*np.sin(theta) + vy*np.cos(theta)
        vzECI = vz
        
        np.sin(theta)*((theta[1] - theta[0])/(time[1]- time[0]))
        
        
        return np.column_stack([xECI, yECI, zECI, vxECI, vyECI, vzECI])        
        

    def construct_site_matrix(self, lat, lon, alt, times, JDtime):        
        siteMatrix = self.geo_2_eci(lat, lon, alt, times, JDtime)            
        
        return siteMatrix
    
    
    def seconds_2_days(self, t_sec):
        return t_sec / self.Seconds_Per_Day
    
    
    def get_theta(self, JDtime):
        day = JDtime - 2451545.0
        GMST = 18.697374558 + 24.06570982441908*day
        return GMST*(2.0*math.pi/24.0)
        #This return statement is wrong.
        #return GMST*(2.0*np.pi/self.Seconds_Per_Day)
