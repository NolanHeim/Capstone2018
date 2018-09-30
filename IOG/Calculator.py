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
from bisect import bisect_left

class Calculator:
    
    def __init__(self):
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
        self.hermiteError = 0.01
        self.timeStepTolerance = 0.01    
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    def cubic_hermite_poly(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti):
        condition = lambda t: tiMinus <= t <= ti        
        Ci = lambda t: ((Vi*(3*hi*((t - tiMinus)**2) - 2*((t - tiMinus)**2)/(hi**3))) +
                        (ViMinus*(hi**3 - 3*hi*((t - tiMinus)**2) + 2*((t - tiMinus)**2))/(hi**3)) +
                        (dVi*(((t - tiMinus)**2)*(t - ti))/(hi**2)) +
                        (dViMinus*((t - tiMinus)*((t - ti)**2))/(hi**2)) )
        return [condition, Ci]
        
    def max_time_step(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a5 = self.compute_a5(hi, ViMinus, Vi, dViMinus, dViHalf, dVi)
        a4 = self.compute_a4(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
        maxDenominator = max([(5*a5*t + a4) for t in [tiMinus, ti]])
        hMax = ((16*self.hermiteError)/(maxDenominator))**(1/4)      
        return hMax
        
    def binary_List_Search(self,dataList, target):
        index = bisect_left(dataList, target)
        if index == 0:
            return 0
        if index == len(dataList):
            return len(dataList)
        beforeIndex = dataList[index - 1]
        afterIndex = dataList[index]
        if afterIndex - target < target - beforeIndex:
            return index-1
        else:
            return index       
        
    #This should combine the instances of each piecewise cubic hermite.
    def cubic_hermite_composite(self, data, position):
        #This is the master method
        
        #Initial Conditions
        VF = self.satellite_visibility(data, position)
        times = [row[1] for row in data]
        
        indexMinus = 0        
        tiMinus = times[indexMinus]
        hiMinus = 100 #Might have to tweak this parameter
        hi = hiMinus
        index = self.binary_List_Search(times, hiMinus)
        ti = times[index]
        
        endOfTime = False
        
        #Loop through the potential values
        while endOfTime == False:
            polySlices = []
            conditionSlices = []
            #Compute Derivatives
            ViMinus = (VF[indexMinus+1] - VF[indexMinus])/hiMinus
            ViHalf = (VF[indexMinus+1] - VF[index-1])/hiMinus
            Vi = (VF[index+1] - VF[index])/hiMinus
            
            #Iterate through the max step (Need to add a maximum iteration reached check)
            k = 1    
            kTolMet = False
            while (k < 100 or (hiMinus > hi))  and kTolMet == False:
                hi = self.max_time_step(hi, ViMinus, ViHalf, Vi, tiMinus, ti)                
                expTol = (numpy.abs(hi - hiMinus))/hiMinus
                if(expTol > self.timeStepTolerance):
                    hiMinus = hi
                else:
                    #Condition is met
                    kTolMet = True
                k = k + 1
            
            hiMinus = hi
            #Interpolate based on the time step hi
            tiMinus = ti
            index = indexMinus
            #Might need to add a 'round down' condition to this for accuracy
            index = self.binary_List_Search(times, tiMinus+hi)
            ti = times[index]
            if(ti >= times[-1]):
                ti = times[-1]
                endOfTime = True
            [condition, Ci] = self.cubic_hermite_poly(tiMinus, ti, ViMinus, Vi)
            polySlices.append(Ci)
            conditionSlices.append(condition)
        
        #Combine the polySlices and conditionSlices to form the piecewise 
             #cubic hermite interpolating polynomial
         
        #When using the output of this, need numpy.asscalar(ANS) to retreive value
        cubicHermitePoly = lambda t: numpy.asscalar(numpy.piecewise(t, conditionSlices, polySlices))
        
        return cubicHermitePoly
        
    def compute_a5(self, hi, ViMinus, Vi, dViMinus, dViHalf, dVi):
        a5 = (24/(hi**5))*(ViMinus - Vi) + (4/(hi**4))*(dViMinus + 4*dViHalf + dVi)        
        return a5

    def compute_a4(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a4 = ( (4/(hi**4))*(ViMinus + 4*ViHalf + Vi) -
                (4/(hi**4))*(dViMinus*(2*tiMinus + 3*ti) + 10*dViHalf*(tiMinus - ti) + dVi*(3*tiMinus - 2*ti)) -
                (24/(hi**5))*(ViMinus*(2*tiMinus + 3*ti) - Vi*(3*tiMinus + 2*ti)) )
        return a4
        
    #UNTESTED
    #Converts the satellite position and time data into sin(theta) vs. seconds  
    #Assumes Position is in Lat/Long/Height
    #Assumes data is the satellite data matrix      
    def satellite_visibility(self, data, position):
        times = [row[1] for row in data]
        positionECEF = self.geodetic_to_ECEF(position[0], position[1], position[2])
        angles = [((math.pi/2) - self.angle_between_vectors(line[1:3], positionECEF)) for line in data]
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
    def angle_between_vectors(self, V1, V2):
        magnitude1 = math.sqrt(V1[0]**2 + V1[1]**2 + V1[2]**2)
        magnitude2 = math.sqrt(V2[0]**2 + V2[1]**2 + V2[2]**2)
        dot = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]
        theta = math.acos(dot/(magnitude1*magnitude2))
        
        return theta
        