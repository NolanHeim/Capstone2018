# -*- coding: utf-8 -*-
#
# CLASSNAME.py
#
# DESCRIPTION
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
#For numpy.piecewise
import numpy
import math

class Calculator:
    
    def __init__(self):
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    def cubicHermitePolynomial(self, subt1, subt2, func):
        condition = lambda t: subt1 <= t <= subt2        
        h = subt2-subt1
        f1 = func[subt1]
        f2 = func[subt2]
        Ci = lambda t: h*f1*f2*t#Function Definition Here
        
        return [condition, Ci]
        
    def visibilityFunction(self, t1, t2, func, N):
        print('Xi')

    #UNTESTED
    #Converts the satellite position and time data into sin(theta) vs. seconds  
    #Assumes Position is in Lat/Long/Height
    #Assumes data is the satellite data matrix      
    def satelliteVisibility(self, data, position):
        times = [row[1] for row in data]
        positionECEF = self.geodetic_to_ECEF(position[0], position[1], position[2])
        angles = [((math.pi/2) - self.angleBetweenVectors(line[1:3], positionECEF)) for line in data]
        sinAngles = math.sin(angles)
        
        return [times, sinAngles]
        
    #Converts geodetic coordinates to ECEF
    def geodetic_to_ECEF(self, Lat, Long, height):
        a = self.equitorialRadius
        b = self.polarRadius
        phi = math.radians(Lat)
        lam = math.radians(Long)
        N = (a**2)/(math.sqrt( ((a**2)*(math.cos(phi)**2)) + ((b**2)*(math.sin(phi)**2)) ))
        X = (N + height)*math.cos(phi)*math.cos(lam)
        Y = (N + height)*math.cos(phi)*math.sin(lam)
        Z = (((b**2)/(a**2))*N + height)*math.sin(phi)
        
        return [X,Y,Z]
        
    #Returns the angle between two vectors in radians
    def angleBetweenVectors(self, V1, V2):
        magnitude1 = math.sqrt(V1[0]**2 + V1[1]**2 + V1[2]**2)
        magnitude2 = math.sqrt(V2[0]**2 + V2[1]**2 + V2[2]**2)
        dot = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]
        theta = math.acos(dot/(magnitude1*magnitude2))
        
        return theta
        
C = Calculator()
a = C.geodetic_to_ECEF(10,10,10)
b = C.angleBetweenVectors(a, [1,1,1])