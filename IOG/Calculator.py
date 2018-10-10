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
import numpy as np
import math
import matplotlib.pyplot as plt
from Transformer import *
from bisect import bisect_left
import time as sleeper
from mpl_toolkits.mplot3d import Axes3D


class Calculator:
    
    def __init__(self):
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
        self.hermiteError = 0.1
        self.timeStepTolerance = 0.05
        self.Rotational_Speed_Earth = (7.2921159 * np.power(10.0, -5.0))
        self.plot = False
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    
    def generate_imaging_opportunities(self, mission, dataMatrices):
        missionCoordinates = mission.get_coordinates()
        #missionStart = mission.get_interval_start_time()
        #missionEnd = mission.get_interval_end_time()
        [mStartJDTime, mEndJDTime] = mission.get_JD_time()
        position = np.array(missionCoordinates[0])
        print(position)
        #Check to see if there is only 1 set of satellite data (matrix is only 2D)
        if(len(dataMatrices.shape) > 2):       
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
        else:
            poly = self.cubic_hermite_composite(dataMatrices, position, mStartJDTime)
            times = dataMatrices[:,0]/3600.0
            y = [poly(t) for t in times]
            plt.plot(times, y)
            plt.show()
            
    def cubic_hermite_poly(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti):
        condition = lambda x: (tiMinus <= x) & (x < ti)        
        poly = lambda x: ( (Vi*(3.0*hi*np.power((x - tiMinus),2.0) - 2.0*np.power((x - tiMinus),2.0)/np.power(hi,3.0))) +
                        (ViMinus*(hi**3.0 - 3.0*hi*np.power((x - tiMinus),2.0) + 2.0*np.power((x - tiMinus),2.0))/np.power(hi,3.0)) +
                        (dVi*(np.power((x - tiMinus),2.0)*(x - ti))/np.power(hi,2.0)) +
                        (dViMinus*((x - tiMinus)*np.power((x - ti),2.0))/np.power(hi,2.0)) )
        return [condition, poly]
    
    
    def max_time_step(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a5 = self.compute_a5(hi, ViMinus, Vi, dViMinus, dViHalf, dVi)
        a4 = self.compute_a4(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
        maxDenominator = max([np.abs(5.0*a5*t + a4) for t in [tiMinus, ti]])
        #print(maxDenominator)
        hMax = np.power(((16.0*self.hermiteError)/(maxDenominator)),0.25) 
            
        return hMax

        
    def binary_List_Search(self, dataList, target):
        index = (np.abs(dataList - target)).argmin()
        return index

            
    #Computes the derivative of the visibility function at a given index
    #TODO Generate a look up table of these (we will be computing these multiple times.)
    def compute_dV(self, index, dataECI, positionECI):
        r_sat = dataECI[index,0:3]        
        v_sat = dataECI[index,3:6]

        r_site = positionECI[index,0:3]        
        v_site = positionECI[index,3:6]        
        m_v_site = np.sqrt(np.sum(v_site*v_site))    
        
        m_r_site = np.sqrt(np.sum(r_site*r_site))
        r_unit_site = r_site/m_r_site
        
        delta_r = r_sat-r_site        
        m_delta_r = np.sqrt(np.sum(delta_r*delta_r))        
        
        d_delta_r = v_sat - v_site
        v_site_unit = v_site/m_v_site

        dV = ( ((np.sum(d_delta_r*r_unit_site) + np.sum(delta_r*v_site_unit))/m_delta_r) - 
                ((np.sum(delta_r*d_delta_r)*np.sum(delta_r*r_unit_site))/np.power(m_delta_r,3.0)))
        
        return dV
        
    #Dot produce of two 3D vectors (i.e. lists)
    def dot_product(self, v1, v2):
        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
        

    def unit_vector(self, v):
        mag = self.get_vec_magnitude(v)
        return [v[0]/mag, v[1]/mag, v[2]/mag]        
        
    #This should combine the instances of each piecewise cubic hermite.
    def cubic_hermite_composite(self, data, position, JDtime):
        #This is the master method
        
        #Initial Conditions
        times = data[:,0]
        unitConversion = Transformer()
        dataECI = unitConversion.ecef_2_eci(data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6],times, JDtime)        
        positionECI = unitConversion.construct_site_matrix(position[0], position[1], position[2], times, JDtime)

        print("Position Matrix Data: X difference is "+str(positionECI[800][0] - positionECI[0][0]))
        print("They are "+str(positionECI[800][0])+" "+str(positionECI[0][0]))
        print("Position Matrix Dimensions are "+str(len(positionECI))+" "+str(len(positionECI[0])))

        VF = self.satellite_visibility(dataECI, times, positionECI)

        if(self.plot):
            plt.plot(times, VF)
            #fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.plot3D(dataECI[:,0], dataECI[:,1], dataECI[:,2], 'g:')
            #ax.scatter3D(positionECI[:,0], positionECI[:,1], positionECI[:,2], 'bo')
            #ax.plot3D(dataECI[:,3], dataECI[:,4], dataECI[:,5])
            #ax.plot3D(positionECI[:,3], positionECI[:,4], positionECI[:,5], 'g')
            plt.show()

        indexMinus = 0
        tiMinus = times[indexMinus]
        hiMinus = 600000 #Might have to tweak this parameter
                
        hi = hiMinus
        index = self.binary_List_Search(times, hiMinus)
        indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))
        
        print("Minus Half Index: "+str(indexMinus)+" "+str(indexHalf)+" "+str(index))        

        ti = times[index]
        
        if(self.plot):
            endOfTime = True        
        else:
            endOfTime = False
    
        polySlices = []
        conditionSlices = []
        #Loop through the potential values
        while endOfTime == False:
            #Iterate through the max step (Need to add a maximum iteration reached check)
            k = 1    
            kTolMet = False
            while k < 100 and kTolMet == False:
                print(str(indexMinus) + ' ' + str(indexHalf) + ' ' + str(index))
                #Compute Derivatives
                ViMinus = VF[indexMinus]
                ViHalf = VF[indexHalf] #THIS NEEDS TO CHANGE
                Vi = VF[index]
                dViMinus = self.compute_dV(indexMinus, dataECI, positionECI)
                dViHalf = self.compute_dV(indexHalf, dataECI, positionECI)
                dVi = self.compute_dV(index, dataECI, positionECI)
                #print([ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi])                
                
                
                hi = self.max_time_step(hiMinus, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
                
                #print('Max TIME')                
                
                sleeper.sleep(0.01)
                expTol = (np.abs(hi - hiMinus))/hiMinus

                if(expTol > self.timeStepTolerance and hi > hiMinus):
                    hiMinus = hi
                elif (hi < hiMinus):
                    kTolMet = True
                    tiMinus = ti
                    indexMinus = index
                else:
                    #Condition is met
                    kTolMet = True
                    hiMinus = hi
                    tiMinus = ti
                    indexMinus = index
            
                index = self.binary_List_Search(times, tiMinus+hi)
                indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))
                ti = times[index]
                print(str(indexMinus) + ' ' + str(indexHalf) + ' ' + str(index))
                k = k + 1
            print('----Done')
            #Need to fix this, right now the start is overridden
            sleeper.sleep(0.05)
            print(hi)
            #Interpolate based on the time step hi
            #Might need to add a 'round down' condition to this for accuracy
            if(ti >= times[-1]):
                ti = times[-1]
                endOfTime = True
            [condition, poly] = self.cubic_hermite_poly(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti)
            conditionSlices.append(condition)            
            polySlices.append(poly)
            
        
        #Combine the polySlices and conditionSlices to form the piecewise 
             #cubic hermite interpolating polynomial
         
        cubicHermitePoly = lambda x: self.piecewise(x, conditionSlices, polySlices)
        #cubicHermitePoly = np.select(conditionSlices(times), polySlices(times))
        print(cubicHermitePoly([times[5000]]))
        return cubicHermitePoly
        
    def compute_a5(self, hi, ViMinus, Vi, dViMinus, dViHalf, dVi):
        a5 = (24.0/(np.power(hi,5.0)))*(ViMinus - Vi) + (4.0/(np.power(hi,4.0)))*(dViMinus + 4.0*dViHalf + dVi)        
        return a5

    def compute_a4(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a4 = ( (4.0/(np.power(hi,4.0)))*(ViMinus + 4.0*ViHalf + Vi) -
                (4.0/(np.power(hi,4.0)))*(dViMinus*(2.0*tiMinus + 3.0*ti) + 10.0*dViHalf*(tiMinus + ti) + dVi*(3.0*tiMinus + 2.0*ti)) -
                (24.0/(np.power(hi,5.0)))*(ViMinus*(2.0*tiMinus + 3.0*ti) - Vi*(3.0*tiMinus + 2.0*ti)) )
        return a4

        
    #UNTESTED
    #Converts the satellite position and time data into sin(theta) vs. seconds  
    #Assumes data is the satellite data matrix
    #Assumes position is the site position matrix
    def satellite_visibility(self, dataECI, times, positionECI):  
        r_sat = dataECI[:,0:3]        

        r_site = positionECI[:,0:3]
        delta_r = r_sat - r_site
        m_delta_r = np.sqrt(np.sum(delta_r*delta_r, axis=1))
        
        print("Size of mag of delta r is "+str(len(m_delta_r)))        

        m_r_site = np.sqrt(np.sum(r_site*r_site, axis=1))
        r_unit_site = r_site/m_r_site[:,None]

        print("Size of r unit site is "+str(len(r_unit_site)))        
        print("Size of r unit site [0] is "+str(len(r_unit_site[0])))
        print(str(r_unit_site[0]))

        numerator = np.sum(delta_r*r_unit_site, axis=1)   
        
        VF = numerator/m_delta_r        
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
        
    def get_visibility_cone(self, Lat, Long, dataECI, times, positionECI, theta0):   
        Lat = np.radians(Lat)
        Long = np.radians(Long)        

        r_sat = dataECI[:,0:3]
        v_sat = dataECI[:,3:6]
        
        r_site = positionECI[:,0:3]
        m_r_site = np.sqrt(np.sum(r_site*r_site, axis=1))
        
        p = np.cross(r_sat, v_sat)
        m_q = m_r_site.max()
        m_p = np.sqrt(np.sum(p*p, axis=1))
        p_tilde = p/m_p[:,None]

        gamma = theta0 + np.arcsin((m_r_site*np.sin((np.pi/2.0) + theta0))/(m_q))
        
        outOfTimeWindow = False
        maxTime = times.max()
        tInOutArray = []
        m = 0
        
        #I don't understand this.
        while(outOfTimeWindow == False):
            tInOut = ((1/self.Rotational_Speed_Earth)*(np.arcsin((np.cos(gamma) - 
                    p_tilde[:,2]*np.sin(Lat))/(np.sqrt(p_tilde[:,0]*p_tilde[:,0] + 
                    p_tilde[:,1]*p_tilde[:,1])*np.cos(Lat))) - Long - 
                    np.arctan(p_tilde[:,0]/p_tilde[:,1]) + 2*np.pi*m))
            #if(tInOut.max() > maxTime):
            #    outOfTimeWindow = True
            #else:
            tInOutArray.append(tInOut)
            outOfTimeWindow = True
            #    m = m + 1
        
        return np.array(tInOutArray)
        
    def piecewise(self, x, conditions, functions):
        conditions = np.asarray(conditions)
        functions = np.asarray(functions)
        evalFcn = []
        evalCond = np.array([cond(x) for cond in conditions])
        #To allow for x to be a vector:        
        for xi in range(0,len(x)):
            yi = np.transpose(evalCond)[xi, :]
            print(yi)
            print(functions[yi])
            evalFcn.append(functions[yi][0](0))
      
        return np.array(evalFcn)
        
        