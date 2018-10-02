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
# TODO: Cite all equations from paper & paper itself.
import numpy
import math
import matplotlib.pyplot as plt
from Transformer import *
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
        #missionStart = mission.get_interval_start_time()
        #missionEnd = mission.get_interval_end_time()
        [mStartJDTime, mEndJDTime] = mission.get_JD_time()
        position = missionCoordinates[0]
        position.append(0)
        print(position)
        
        for dataMatrix in dataMatrices:
            poly = self.cubic_hermite_composite(dataMatrix, position, mStartJDTime)
            print(poly)
            t = [row[0] for row in dataMatrix]
            y = [poly[ti] for ti in t]
            plt.plot(t, y)
            plt.show()
            VF = self.satellite_visibility(dataMatrix, position, mStartJDTime)
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
            
    #Computes the derivative of the visibility function at a given index
    def compute_dV(self, index, dataECI, positionECI):
        r_sat = [dataECI[index][0], dataECI[index][1], dataECI[index][2]]
        v_sat = [dataECI[index][3], dataECI[index][4], dataECI[index][5]]
                
        r_site = [positionECI[0], positionECI[1], positionECI[2]]
        r_unit_site = [(i/self.get_vec_magnitude(r_site)) for i in r_site]

        delta_r = [r_sat[0]-r_site[0], r_sat[1]-r_site[1], r_sat[2]-r_site[2]]
        m_delta_r = self.get_vec_magnitude(delta_r)
        
        #What is v_site? zero?
        d_delta_r = v_sat
        v_site_unit = [0,0,0]
        
        dV = (((self.dot_product(d_delta_r, r_unit_site) + self.dot_product(delta_r, v_site_unit))/(m_delta_r)) - 
            ((self.dot_product(delta_r, d_delta_r)*self.dot_product(delta_r, r_unit_site))/(m_delta_r**3.0)))
        
        return dV
        
    #Dot produce of two 3D vectors (i.e. lists)
    def dot_product(self, v1, v2):
        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
        
    #This should combine the instances of each piecewise cubic hermite.
    def cubic_hermite_composite(self, data, position, JDtime):
        #This is the master method
        
        #Initial Conditions
        times = [row[0] for row in data]     
        unitConversion = Transformer()
        dataECI = [unitConversion.ecef_2_eci(l[1],l[2],l[3],l[4],l[5],l[6], JDtime) for l in data]
        #MOVE Computation of position in ECI to here
        positionECI = unitConversion.geo_2_eci(position[0], position[1], position[2], JDtime)
        VF = self.satellite_visibility(dataECI, times, positionECI)
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
            dViMinus = self.compute_dV(indexMinus, dataECI, positionECI)
            dViHalf = self.compute_dV(indexMinus+1, dataECI, positionECI) # TODO: Change Index
            dVi = self.compute_dV(index, dataECI, positionECI)
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
    #Assumes data is the satellite data matrix      
    def satellite_visibility(self, dataECI, times, positionECI):      
        r_sat = [[line[0], line[1], line[2]] for line in dataECI]
        r_site = [positionECI[0], positionECI[1], positionECI[2]]

        delta_r = [[line[0]-r_site[0], line[1]-r_site[1], line[2]-r_site[2]] for line in r_sat]
        m_delta_r = [self.get_vec_magnitude(line) for line in delta_r]

        r_unit_site = [(i/self.get_vec_magnitude(r_site)) for i in r_site]

        numerator = [self.dot_product(line, r_unit_site) for line in delta_r]
        
        VF = []

        for index in range(0,len(numerator)):
            VF.append(numerator[index]/m_delta_r[index])
        
        return VF
    
    
    
    
    #WILL LIKELY BE DEPRECIATED AS ECEF IS NOT THE WORKING COORDINATE SYSTEM
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
    
    #returns the magnitude of the given vector    
    def get_vec_magnitude(self, V):
        magnitude = math.sqrt(V[0]**2 + V[1]**2 + V[2]**2)
        
        return magnitude