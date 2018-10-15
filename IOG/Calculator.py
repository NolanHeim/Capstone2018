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
import matplotlib.pyplot as plt
import time as sleeper
from mpl_toolkits.mplot3d import Axes3D


class Calculator:
    
    def __init__(self):
        self.hermiteError = 0.01
        self.timeStepTolerance = 0.05
        self.plot = False #To display the resulting plots
        self.verbose = False #To display extra information at each step.
        self.maxIterations = 100
        self.initialTimeStep = 120
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
    
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    def generate_imaging_opportunities(self, mission, dataMatrices):
        missionCoordinates = mission.get_coordinates()
        position = np.array(missionCoordinates[0])
        if(self.verbose):        
            print(position)
        
        #Check to see if there is only 1 set of satellite data (matrix is only 2D)
        if(len(dataMatrices.shape) > 2):       
            for dataMatrix in dataMatrices:
                [poly, timeWindows] = self.cubic_hermite_composite(dataMatrix, position)
                times = dataMatrix[:,0]
                y = poly(times)
                if(self.plot):
                    plt.plot(times, y)
                    plt.show()

        else:
            [poly, timeWindows] = self.cubic_hermite_composite(dataMatrices, position)
            times = dataMatrices[:,0]
            if(self.plot):
                dataECEF = dataMatrices[:,1:7]
                posECEF = self.geo_to_ecef(position[0], position[1], position[2])
                VF = self.satellite_visibility(dataECEF, times, posECEF)                        
                y = poly(times)
                fig1, ax1 = plt.subplots()
                ax1.plot(times, y)
                ax1.plot(times, VF)
                print(timeWindows)
                plt.show()
                
        
    def cubic_hermite_poly(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti):
        condition = lambda x: (tiMinus <= x) & (x < ti)
        if(self.verbose):        
            print([tiMinus,ti])
        poly = lambda x: ( (Vi*((3.0*hi*np.power((x - tiMinus),2.0) - 2.0*np.power((x - tiMinus),3.0))/np.power(hi,3.0))) +
                        (ViMinus*(np.power(hi,3.0) - 3.0*hi*np.power((x - tiMinus),2.0) + 2.0*np.power((x - tiMinus),3.0))/np.power(hi,3.0)) +
                        (dVi*(np.power((x - tiMinus),2.0)*(x - ti))/np.power(hi,2.0)) +
                        (dViMinus*((x - tiMinus)*np.power((x - ti),2.0))/np.power(hi,2.0)) )
        
        roots = self.compute_cubic_roots(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti, condition)
        return [condition, poly, roots]
    
    
    def compute_cubic_roots(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti, domain):
        #Third Order Coefficient
        t3 = (1.0/np.power(hi, 3.0))*(-2.0*Vi + 2.0*ViMinus + hi*dVi + hi*dViMinus)
        
        t2 = (1.0/np.power(hi, 3.0))*( (3.0*hi + 6.0*tiMinus)*Vi 
                                    - (3.0*hi + 6.0*tiMinus)*ViMinus 
                                    - (hi*(ti + 2.0*tiMinus))*dVi
                                    - (hi*(2.0*ti + tiMinus))*dViMinus )
        
        t1 = (1.0/np.power(hi, 3.0))*( -1.0*(6.0*hi*tiMinus + 6.0*np.power(tiMinus, 2.0))*Vi 
                                    + (6.0*hi*tiMinus + 6.0*np.power(tiMinus, 2.0))*ViMinus 
                                    + hi*tiMinus*(2.0*ti + tiMinus)*dVi
                                    + hi*ti*(ti + 2.0*tiMinus)*dViMinus )
        
        t0 = (1.0/np.power(hi, 3.0))*( (3.0*hi*np.power(tiMinus, 2.0) + 2.0*np.power(tiMinus, 3.0))*Vi
                + (np.power(hi, 3.0) - 3.0*hi*np.power(tiMinus, 2.0) - 2.0*np.power(tiMinus, 3.0))*ViMinus
                - hi*ti*np.power(tiMinus, 2.0)*dVi
                - hi*tiMinus*np.power(ti, 2.0)*dViMinus )
        
        coefficients = np.array([t3,t2,t1,t0])
        complexRoots = np.roots(coefficients)
        #Take only real roots
        realRoots = np.isreal(complexRoots)
        roots = np.real(complexRoots[realRoots])
        
        inDomain = domain(roots)
        
        if(self.verbose):        
            print([tiMinus, roots, ti])                
            print(inDomain)
            sleeper.sleep(0.2)

        validRoots = roots[inDomain]
        
        return validRoots
    
    
    def max_time_step(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a5 = self.compute_a5(hi, ViMinus, Vi, dViMinus, dViHalf, dVi)
        a4 = self.compute_a4(hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
        maxDenominator = max([np.abs(5.0*a5*t + a4) for t in [tiMinus, ti]])

        hMax = np.power(((16.0*self.hermiteError)/(maxDenominator)),0.25) 
            
        return hMax

        
    def binary_List_Search(self, dataList, target):
        index = (np.abs(dataList - target)).argmin()
        return index


    def compute_dV(self, dataECEF, posECEF):
        r_sat = dataECEF[:,0:3]        
        v_sat = dataECEF[:,3:6]

        r_site = posECEF      
        
        m_r_site = np.sqrt(np.sum(r_site*r_site, axis=1))
        r_unit_site = r_site/m_r_site[:,None]
        
        delta_r = r_sat-r_site        
        m_delta_r = np.sqrt(np.sum(delta_r*delta_r, axis=1))

        d_delta_r = v_sat
        
        dV = (np.sum(d_delta_r*r_unit_site, axis=1))/m_delta_r
     
        term2 = ((np.sum(delta_r*d_delta_r, axis=1)*np.sum(delta_r*r_unit_site, axis=1)))

        term2Denominator = np.power(m_delta_r,-3.0)
        term2 *=term2Denominator
        dV -= term2
        
        return dV
    
    
    #This will combine the instances of each piecewise cubic hermite.
    def cubic_hermite_composite(self, data, position):        
        #Initial Conditions
        times = data[:,0]
        dataECEF = data[:,1:7]
        posECEF = self.geo_to_ecef(position[0], position[1], position[2])

        VF = self.satellite_visibility(dataECEF, times, posECEF)
        dVF = self.compute_dV(dataECEF, posECEF)
        
        if(self.plot):
            plt.plot(times, VF)
            #fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.plot3D(dataECEF[:,0], dataECEF[:,1], dataECEF[:,2], 'g:') #Position
            #ax.scatter3D(posECEF[:,0], posECEF[:,1], posECEF[:,2], 'bo') #Position
            #ax.plot3D(dataECEF[:,3], dataECEF[:,4], dataECEF[:,5]) #Velocity
            plt.show()

        indexMinus = 0
        tiMinus = times[indexMinus]
        hiMinus = self.initialTimeStep #Might have to tweak this parameter
                
        hi = hiMinus
        index = self.binary_List_Search(times, hiMinus)
        indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))
        
        if(self.verbose):
            print("MinusIndex, HalfIndex, Index: "+str(indexMinus)+", "+str(indexHalf)+", "+str(index))        

        ti = times[index]
    
        polySlices = []
        conditionSlices = []
        rootSlices = []
        endOfTime= False        
        #Loop through the potential values
        while endOfTime == False:
            #Iterate through the max step (Need to add a maximum iteration reached check)
            k = 1    
            kTolMet = False
            while k < self.maxIterations and kTolMet == False:
                #Compute Visibility Function and Derivatives.
                ViMinus = VF[indexMinus]
                ViHalf = VF[indexHalf]
                Vi = VF[index]
                dViMinus = dVF[indexMinus]
                dViHalf = dVF[indexHalf]
                dVi = dVF[index]           
                
                hi = self.max_time_step(hiMinus, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti)
                expTol = (np.abs(hi - hiMinus))/hiMinus

                if(expTol > self.timeStepTolerance and hi > hiMinus):
                    hiMinus = hi
                    index = self.binary_List_Search(times, tiMinus+hi)
                    indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))
                    ti = times[index]
                elif(hi < hiMinus):
                    hi = hiMinus
                    kTolMet = True
                else:
                    #Condition is met
                    kTolMet = True
                    
                k = k + 1
                if(self.verbose):
                    print(str(indexMinus) + ' ' + str(indexHalf) + ' ' + str(index))
            
                
            if(ti >= times[-1]):
                ti = times[-1]
                endOfTime = True
            
            [condition, poly, roots] = self.cubic_hermite_poly(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti)
            conditionSlices.append(condition)            
            polySlices.append(poly)            
            rootSlices.append(roots)
            
            if(self.verbose):
                print('hi: ' + str(hi))
                sleeper.sleep(2)
            
            hiMinus = hi
            tiMinus = ti
            indexMinus = index            
            index = self.binary_List_Search(times, tiMinus+hi)
            indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))     
            ti = times[index]

        #Piecewise cubic hermite interpolating polynomial
        timeWindows = self.create_time_windows(rootSlices, times, VF)
        cubicHermitePoly = lambda x: self.piecewise(x, conditionSlices, polySlices)
        #print(cubicHermitePoly([100,10000]))
        return [cubicHermitePoly, timeWindows]    
    
    
    def create_time_windows(self, roots, times, VF):
        
        rootsList = []        
        #Process the roots list
        for interval in roots:
            for root in interval:
                rootsList.append(root)
        
        rootsList = np.array(rootsList)
        if(self.verbose):        
            print(rootsList)
        timingWindow = []
        startIndex = 0
        
        if(VF[0] > 0):  #then we are already in visibility range, so t0 to the first root is an interval
            timingWindow.append([0, rootsList[0]])
            startIndex = 1
        elif(VF[0] < 0):  #then we are not in visibility range, so the first root will be the start of an interval
            timingWindow.append([rootsList[0], rootsList[1]])
            startIndex = 2
        else:            
            #It is zero at t = 0, check if star/end of window
            if(rootsList[0] == 0.0):
                indexHalf = self.binary_List_Search(times, rootsList[1]/2.0)
                endrootIndex = 1
            else:
                indexHalf = self.binary_List_Search(times, rootsList[0]/2.0)
                endrootIndex = 0
            
            if(VF[indexHalf] > 0):
                timingWindow = [0, rootsList[endrootIndex]]
                startIndex = endrootIndex+1
            else:
                timingWindow = [rootsList[endrootIndex], rootsList[endrootIndex+1]]
                startIndex = endrootIndex+2
        
#        endIndex = len(rootsList)+1
        endIndex = len(rootsList)
        for index in range(startIndex, endIndex, 2):
            timingWindow.append([rootsList[index], rootsList[index+1]])
        
        #Check to see if the roots form a closed set.
        
        if(range(startIndex, endIndex)[-1] == (endIndex-2)):          
            if(rootsList[-1] != times[-1]):
                timingWindow.append([rootsList[-1], times[-1]])
        
        if(self.verbose):
            print(timingWindow)
        
        return np.array(timingWindow)
    
    
    def compute_a5(self, hi, ViMinus, Vi, dViMinus, dViHalf, dVi):
        a5 = (24.0/(np.power(hi,5.0)))*(ViMinus - Vi) + (4.0/(np.power(hi,4.0)))*(dViMinus + 4.0*dViHalf + dVi)        
        return a5

    
    def compute_a4(self, hi, ViMinus, ViHalf, Vi, dViMinus, dViHalf, dVi, tiMinus, ti):
        a4 = ( (4.0/(np.power(hi,4.0)))*(ViMinus + 4.0*ViHalf + Vi) -
                (4.0/(np.power(hi,4.0)))*(dViMinus*(2.0*tiMinus + 3.0*ti) + 10.0*dViHalf*(tiMinus + ti) + dVi*(3.0*tiMinus + 2.0*ti)) -
                (24.0/(np.power(hi,5.0)))*(ViMinus*(2.0*tiMinus + 3.0*ti) - Vi*(3.0*tiMinus + 2.0*ti)) )
        return a4

        
    def satellite_visibility(self, dataECI, times, positionECI):  
        r_sat = dataECI[:,0:3]        
        r_site = positionECI
        
        delta_r = r_sat - r_site
        m_delta_r = np.sqrt(np.sum(delta_r*delta_r, axis=1))
        
        m_r_site = np.sqrt(np.sum(r_site*r_site, axis=1))
        r_unit_site = r_site/m_r_site[:,None]

        numerator = np.sum(delta_r*r_unit_site, axis=1)   
        
        VF = numerator/m_delta_r        
        return VF


    def piecewise(self, x, conditions, functions):
        conditions = np.asarray(conditions)
        functions = np.asarray(functions)
        evalFcn = []
        evalCond = np.array([cond(x) for cond in conditions])
        #To allow for x to be a vector:        
        for xi in range(0,(len(x)-1)):
            yi = np.transpose(evalCond)[xi, :]
            evalFcn.append(functions[yi][0](x[xi]))
        #Add the final point
        yf = np.transpose(evalCond)[len(x)-2, :]
        evalFcn.append(functions[yf][0](x[-1]))
        
        return np.array(evalFcn)


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