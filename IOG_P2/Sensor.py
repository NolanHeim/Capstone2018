# -*- coding: utf-8 -*-
#
# Sensor.py
#
# Sensor Calculations
#
# Initial Creation Date: 11/17/2018
#
# Written by Jordan Jones and Nolan Heim
#

import numpy as np
import math
from Calculator import *

class Sensor:
    
    def __init__(self):
        self.start = True
        self.calculator = Calculator()
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
    
    ####Need consistent satview rectangle points
    ###[[lambda1, phi1],[lambda1,phi2],[lambda2,phi1],[lambda2,phi2]]
    #Determines the intersection between the satellite and
    #the N-D polygon (assume arguments in geo & are numpy arrays)
    def determineIntersection(self, satView, Npoly, satVel):
        #Determine if all points are at one time within
        #the viewing platform
        
        
        #Assuming Rectangle
        #Centriod of satelliteView
        #centriods = (satView[:,0,:] + satView[:,3,:])/2.0
            
            
        #Construct Fn,ij -> Store in dictionaries
        Fnij = {};
        for npoint in range(0, Npoly.shape[0]):
            vertex = Npoly[npoint]
            
            ##NOW Need a function to determine the roots of these functions.
            #Solution: Use the Cubic Hermite Polynomial
            #Get the satellite path
            #Trim the satellite data matrix
            #Send information into 'cubic_hermite_composite'
            #Get back the timing windows & polynomial handle
            
            #Trim the datamatrix for the satellite to only contain the 
            #timing windows.            
            
            #This code needs the trimmed data matrixes for satView.
            D01 = self.getMagnitude(satView[0] - satView[1])
            D02 = self.getMagnitude(satView[0] - satView[2])
            D13 = self.getMagnitude(satView[1] - satView[3])
            D23 = self.getMagnitude(satView[2] - satView[3])
            
            dn0 = self.getMagnitude(vertex - satView[0])
            dn1 = self.getMagnitude(vertex - satView[0])
            dn2 = self.getMagnitude(vertex - satView[0])
            dn3 = self.getMagnitude(vertex - satView[0])
            
            Fnij[npoint] = {}
            Fnij[npoint]['01'] = dn0 + dn1 - D01 
            Fnij[npoint]['02'] = dn0 + dn2 - D02
            Fnij[npoint]['13'] = dn1 + dn3 - D13
            Fnij[npoint]['23'] = dn2 + dn3 - D23            
            
            #Now need to apply the cubic hermite to each of the timing windows.
            #We shall have to see how this goes to determine if we keep use the C.H.P
            #for this part.
            
            #Basically take cubic_hermite_composite(timingData, vertex, epoch)
            #Get back a polynomial for the vertex & it's timing window.
            
        #Stitch all of the timing windows together.
            
            
            
            
    def geo_to_ecef(self, lat, lon, alt):
        a = self.equitorialRadius
        b = self.polarRadius

        phi = np.radians(lat)
        lam = np.radians(lon)
        N = ((np.power(a,2.0))/
            (np.sqrt( np.power(a*np.cos(phi),2.0) + np.power(b*np.sin(phi),2.0))) )
        
        X = (N + alt)*np.cos(phi)*np.cos(lam)
        Y = (N + alt)*np.cos(phi)*np.sin(lam)
        Z = (( np.power(b,2.0)/np.power(a,2.0) )*N + alt)*np.sin(phi)
        
        return np.column_stack([X,Y,Z])
        
     
    #based off the second response here https://gis.stackexchange.com/questions/265909/converting-from-ecef-to-geodetic-coordinates
    def ecef_to_geo(self, position):
        a = self.equitorialRadius
        b = self.polarRadius
        
        lam = np.atan2(y,x)

        e = np.sqrt(np.power(a,2.0)-np.power(b,2.0))/a
        p = np.sqrt(np.power(x,2.0)+np.power(y,2.0))
        
        # have to iterate to find alt, with an initial guess of 0
        alt_0 = 0.0
        phi = math.atan2(z,p*(1.0-np.power(e,2.0)))
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        N = np.power(a,2.0)/np.sqrt(np.power(a*cos_phi,2.0)+np.power(b*sin_phi,2.0))
        alt = p/cos_phi - N
        while abs(alt-alt_0) > 1.0e-6:
            alt_0 = alt
            phi = np.atan2(z,p*(1.0-np.power(e,2.0)*N/(N+alt)))
            cos_phi = np.cos(phi)
            sin_phi = math.sin(phi)
            N = np.power(a,2.0)/np.sqrt(np.power(a*cos_phi,2.0)+np.power(b*sin_phi,2.0))
            alt = p/cos_phi - N

        lon = np.degrees(lam)
        lat = np.degrees(phi)

        geo = np.column_stack([lat, lon, alt])
        
        return geo

    
    def getMagnitude(self, vec):
        return np.sqrt(np.sum(vec, axis=1)) 
        
        
    def viewing_rectangle(self, satelliteData, sensor):
        
        


    def rectangle_center(self, satelliteData, sensor):
        #assume these are in radians        
        phi_h = sensor.get("phi_h")
        phi_v = sensor.get("phi_v")
        d_phi_h = sensor.get("d_phi_h")
        d_phi_v = sensor.get("d_phi_v")
        
        positions = satelliteData[:][1:4]
        velocities = satelliteData[:][4:7]
        
        #these are in degrees (or degrees/sec)
        [lat, lon, alt] = self.ecef_to_geo(positions)
        [v_lat, v_lon, v_alt] = self.ecef_to_geo(velocities)
        
        #in radians
        psi_h = self.get_psi(phi_h, alt)
        psi_v = self.get_psi(phi_v, alt)
        d_psi_h = self.get_psi(d_phi_h, alt)
        d_psi_v = self.get_psi(d_phi_v, alt)
        
        #angle between psi space and actual space (in radians)
        xi = np.atan2(v_lat/v_lon)
        
        #rotate the center point of the viewing rectangle by -xi
        tau = np.cos(xi)*psi_h + np.sin(xi)*psi_v
        omega = -np.sin(xi)*psi_h + np.cos(xi)*psi_v
        
        d_lat = 
        #construct four points
        #call them tau (psi_h) and omega (psi_v)         
        #tau1 = psi_h - d_psi_h
        #tau2 = psi_h + d_psi_h
        #omega1 = psi_v - d_psi_v
        #omega2 = psi_v + d_psi_v
        
        d_lat_rad_1 = np.cos(xi)*tau_1 +   
        d_lat_rad_2 = np.cos(tau2)
        d_lon_rad_1 = 
        d_lon_rad_2 = 
        
                
    def get_psi(self, phi, altitude):
        theta = np.asin((altitude + self.r_earth)*(np.sin(phi)/self.r_earth))
        return (np.pi - theta - phi)
        

        
        
    def viewing_circle():