# -*- coding: utf-8 -*-
#
# CLASSNAME.py
#
# DESCRIPTION
#
# Initial Creation Date: 09/26/2018
# REF: http://aa.usno.navy.mil/faq/docs/GAST.php
# Written by Jordan Jones and Nolan Heim
#
#For numpy.piecewise
import numpy
import math
import matplotlib.pyplot as plt
from bisect import bisect_left
import time as sleeper

class Calculator:
    
    def __init__(self):
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
        self.hermiteError = 0.1
        self.timeStepTolerance = 0.01    
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    
    def generate_imaging_opportunities(self, mission, dataMatrices):
        missionCoordinates = mission.get_coordinates()
        missionStart = mission.get_interval_start_time()
        missionEnd = mission.get_interval_end_time()
        position = missionCoordinates[0]
        position.append(0)
        print(position)
        
        for dataMatrix in dataMatrices:
            poly = self.cubic_hermite_composite(dataMatrix, position)
            print(poly)
            t = [row[0] for row in dataMatrix]
            y = [poly[ti] for ti in t]
            plt.plot(t, y)
            plt.show()
            VF = self.satellite_visibility(dataMatrix, position)
            plt.plot(VF)
            plt.show()
            
            
    
    
    
    def cubic_hermite_poly(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti):
        condition = lambda t: tiMinus <= t <= ti        
        Ci = lambda t: ((Vi*(3.0*hi*((t - tiMinus)**2.0) - 2.0*((t - tiMinus)**2.0)/(hi**3.0))) +
                        (ViMinus*(hi**3.0 - 3.0*hi*((t - tiMinus)**2.0) + 2.0*((t - tiMinus)**2.0))/(hi**3.0)) +
                        (dVi*(((t - tiMinus)**2.0)*(t - ti))/(hi**2.0)) +
                        (dViMinus*((t - tiMinus)*((t - ti)**2))/(hi**2.0)) )
        return [condition, Ci]
        
    def max_time_step(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a5 = self.compute_a5(hi, ViMinus, Vi, dViMinus, dViHalf, dVi)
        a4 = self.compute_a4(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
        print(a4)
        print(a5)
        maxDenominator = max([numpy.abs(5.0*a5*t + a4) for t in [tiMinus, ti]])
        print(maxDenominator)
        hMax = ((16*self.hermiteError)/(maxDenominator))**(0.25) 
            
        return hMax
        
    def binary_List_Search(self,dataList, target):
        index = bisect_left(dataList, target)
        if index == 0:
            return 0
        if index == len(dataList):
            return len(dataList)-1
        beforeIndex = dataList[index]
        afterIndex = dataList[index+1]
        if numpy.abs(afterIndex - target) < numpy.abs(target - beforeIndex):
            return index+1
        else:
            return index     
        
    #This should combine the instances of each piecewise cubic hermite.
    def cubic_hermite_composite(self, data, position):
        #This is the master method
        
        #Initial Conditions
        times = [row[0] for row in data]
        VF = self.satellite_visibility(data, times, position)
        #plt.plot(times, VF)
        #plt.show()
        
        indexMinus = 0        
        tiMinus = times[indexMinus]
        hiMinus = 1200 #Might have to tweak this parameter
        hi = hiMinus
        index = self.binary_List_Search(times, hiMinus)
        ti = times[index]
        
        endOfTime = False
        polySlices = []
        conditionSlices = []
        #Loop through the potential values
        while endOfTime == False:
            
            #Compute Derivatives
            ViMinus = VF[indexMinus]
            ViHalf = (VF[indexMinus] + VF[index])/2 #THIS NEEDS TO CHANGE
            Vi = VF[index]
            dViMinus = (VF[indexMinus+1] - VF[indexMinus])/(times[indexMinus+1] - times[indexMinus])
            dViHalf = (VF[indexMinus] - VF[index])/hiMinus #MIGHT NEED TO CHANGE
            dVi = (VF[index+1] - VF[index])/(times[index+1] - times[index])
            print([ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi])
            #Iterate through the max step (Need to add a maximum iteration reached check)
            k = 1    
            kTolMet = False
            while k < 100 and kTolMet == False:
                hi = self.max_time_step(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti) 
                print('Max TIME')                
                print(hi)
                sleeper.sleep(3)
                expTol = (numpy.abs(hi - hiMinus))/hiMinus
                print('ExpTol')
                print(expTol)
                if(expTol > self.timeStepTolerance):
                    hiMinus = hi
                #elif(hiMinus > hi):
                    #Condition is met
                #    kTolMet = True
                #    hi = hiMinus
                else:
                    #Condition is met
                    kTolMet = True
                k = k + 1
                
            sleeper.sleep(5)
            print(hi)
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
            [condition, Ci] = self.cubic_hermite_poly(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti)
            polySlices.append(Ci)
            conditionSlices.append(condition)
        
        #Combine the polySlices and conditionSlices to form the piecewise 
             #cubic hermite interpolating polynomial
         
        #When using the output of this, need numpy.asscalar(ANS) to retreive value
        cubicHermitePoly = lambda t: numpy.asscalar(numpy.piecewise(t, conditionSlices, polySlices))
        
        return cubicHermitePoly
        
    def compute_a5(self, hi, ViMinus, Vi, dViMinus, dViHalf, dVi):
        a5 = (24.0/(hi**5.0))*(ViMinus - Vi) + (4.0/(hi**4.0))*(dViMinus + 4.0*dViHalf + dVi)        
        return a5

    def compute_a4(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a4 = ( (4.0/(hi**4.0))*(ViMinus + 4.0*ViHalf + Vi) -
                (4.0/(hi**4.0))*(dViMinus*(2.0*tiMinus + 3.0*ti) + 10.0*dViHalf*(tiMinus - ti) + dVi*(3.0*tiMinus - 2.0*ti)) -
                (24.0/(hi**5.0))*(ViMinus*(2.0*tiMinus + 3.0*ti) - Vi*(3.0*tiMinus + 2.0*ti)) )
        return a4
        
    #UNTESTED
    #Converts the satellite position and time data into sin(theta) vs. seconds  
    #Assumes Position is in Lat/Long/Height
    #Assumes data is the satellite data matrix      
    def satellite_visibility(self, data, times, position):
        positionECEF = self.geodetic_to_ECEF(position[0], position[1], position[2])
        angles = [((math.pi/2) - self.angle_between_vectors(line[1:4], positionECEF)) for line in data]     
        sinAngles = [math.sin(angle) for angle in angles]
        
        return sinAngles
        
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
        