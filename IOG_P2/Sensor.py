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
from Calculator import *

class Sensor:
    
    def __init__(self):
        self.start = True
        self.calculator = Calculator()
    
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
            
            
            
            

        
        
    def ecef_to_geo(self, position, velocity):
        return 0
    
    def getMagnitude(self, vec):
        return np.sqrt(np.sum(vec, axis=1)) 
        