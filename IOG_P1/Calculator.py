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
import time
from mpl_toolkits.mplot3d import Axes3D
import datetime


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
        self.Rotational_Speed_Earth = (7.2921159 * np.power(10.0, -5.0))
        self.RataDieJDtime = 1721424.5 #Days
        self.epoch_index = 0
        self.uuid_index = 1
    
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    def generate_imaging_opportunities(self, mission, dataMatrices, extraInfoMatrix):
        #self.t0 = time.time()
        missionCoordinates = mission.get_coordinates()
        position = np.array(missionCoordinates)
        
        mission_interval_start = mission.get_interval_start_time()
        mission_interval_end = mission.get_interval_end_time()
                
        if(self.verbose):        
            print(position)
        
        timingWindows_Matrix = []
        satellites_list = []
        
        for i in range(0, len(dataMatrices)):
            #in the first case, every satellite will be considered
            #in the second case, a given satellite will only be considered if 
            if(mission.get_ids_to_consider() == [] or extraInfoMatrix[i][self.uuid_index] in mission.get_ids_to_consider()):              
                epoch = self.calc_epoch(extraInfoMatrix[i][self.epoch_index])
                times = dataMatrices[i][:,0]            
                
                data_interval_start = epoch
                data_interval_end = epoch + datetime.timedelta(seconds=times[-1])
    
    
                [start_index, end_index] = self.check_time_intersection(mission_interval_start, mission_interval_end, 
                                            data_interval_start, data_interval_end, times)
    
                if(self.verbose):
                    print(mission_interval_start)
                    print(mission_interval_end)
                    print(data_interval_start)
                    print(data_interval_end)            
                    print("Start index: "+str(start_index))
                    print("End index: "+str(end_index))
                
                trimmedMatrix = dataMatrices[i][start_index:end_index]
                
                if(self.verbose):            
                    print(trimmedMatrix.shape)
    
                [poly, timingWindows] = self.cubic_hermite_composite(trimmedMatrix, position, epoch)
    
                #might want to redfinie times here
    
                if(self.plot):
                    stopPlot = 3*1440
                    dataECEF = dataMatrices[i][:,1:7]
                    posECEF = self.geo_to_ecef(position[0], position[1], position[2])
                    VF = self.satellite_visibility(dataECEF, times, posECEF)                        
                    y = poly(times[0:stopPlot])
                    print(len(times))
                    fig1, ax1 = plt.subplots()
                    ax1.plot(times[0:stopPlot], VF[0:stopPlot], color="blue", label="Visibility Function", linewidth=5.0)
                    ax1.plot(times[0:stopPlot], y, color="black", linestyle='-', label="Interpolation", linewidth=3.0)            
                    ax1.set_xlabel('Time (s)')
                    ax1.set_ylabel('Visibility Function')
                    ax1.set_title('Interpolation With Error = 0.1')
                    plt.legend()
                    plt.show()
                
                timingWindows_Matrix.append(timingWindows)
                satellites_list.append(str(extraInfoMatrix[i][self.uuid_index]))

        return [timingWindows_Matrix, satellites_list]
        
        
    def cubic_hermite_poly(self, hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti):
        #t2 = time.time()        
        
        condition = lambda x: (tiMinus <= x) & (x < ti)
        if(self.verbose):        
            print([tiMinus,ti])
            
        poly = lambda x: ( (Vi*((3.0*hi*np.power((x - tiMinus),2.0) - 2.0*np.power((x - tiMinus),3.0))/np.power(hi,3.0))) +
                        (ViMinus*(np.power(hi,3.0) - 3.0*hi*np.power((x - tiMinus),2.0) + 2.0*np.power((x - tiMinus),3.0))/np.power(hi,3.0)) +
                        (dVi*(np.power((x - tiMinus),2.0)*(x - ti))/np.power(hi,2.0)) +
                        (dViMinus*((x - tiMinus)*np.power((x - ti),2.0))/np.power(hi,2.0)) )
        
        roots = np.sort(self.compute_cubic_roots(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti, condition))

        #print("Time for interpolation + root finding = "+str(time.time() - t2))

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
            time.sleep(0.2)

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
    def cubic_hermite_composite(self, data, position, epoch):        
        #self.t1 = time.time()
        #print("Time for data trimming = "+str(self.t1 - self.t0))        
        
        #Initial Conditions
        times = data[:,0]
        dataECEF = data[:,1:7]

        posECEF = self.geo_to_ecef(position[0], position[1], position[2])

        VF = self.satellite_visibility(dataECEF, times, posECEF)
        dVF = self.compute_dV(dataECEF, posECEF)
        
        #if(self.plot):
            #plt.plot(times, VF)
            #fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.plot3D(dataECEF[:,0], dataECEF[:,1], dataECEF[:,2], 'g:') #Position
            #ax.scatter3D(posECEF[:,0], posECEF[:,1], posECEF[:,2], 'bo') #Position
            #ax.plot3D(dataECEF[:,3], dataECEF[:,4], dataECEF[:,5]) #Velocity
            #plt.show()

        indexMinus = 0
        tiMinus = times[indexMinus]
        hiMinus = self.initialTimeStep #Might have to tweak this parameter
                
        hi = hiMinus
        index = self.binary_List_Search(times, hiMinus)
        indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))
        
        if(self.verbose):
            print("MinusIndex, HalfIndex, Index: "+str(indexMinus)+", "+str(indexHalf)+", "+str(index))        

        ti = times[index]

    #Related to viewing cone, untested
        #test = self.satellite_viewing_cone(dataECEF,posECEF,position)
        #print(test)
                
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
            
                #print("Time for determining this timestep= "+str(time.time() - self.t1))
                
            if(ti >= times[-1]):
                ti = times[-1]
                endOfTime = True
            
            [condition, poly, roots] = self.cubic_hermite_poly(hi, ViMinus, Vi, dViMinus, dVi, tiMinus, ti)
            conditionSlices.append(condition)            
            polySlices.append(poly)            
            rootSlices.append(roots)
            
            if(self.verbose):
                print('hi: ' + str(hi))
                time.sleep(2)
            
            hiMinus = hi
            tiMinus = ti
            indexMinus = index            
            index = self.binary_List_Search(times, tiMinus+hi)
            indexHalf = self.binary_List_Search(times, tiMinus+(hi/2.0))     
            ti = times[index]
        
        if(self.verbose):        
            print("Poly slices has length: ")
            print(len(polySlices))
       
        #Piecewise cubic hermite interpolating polynomial
        timeWindows = self.create_time_windows(rootSlices, times, VF, epoch)

        if(self.verbose):        
            print(len(timeWindows))
        
        cubicHermitePoly = lambda x: self.piecewise(x, conditionSlices, polySlices)
        #print(cubicHermitePoly([100,10000]))
        return [cubicHermitePoly, timeWindows]    
    
    
    def create_time_windows(self, roots, times, VF, epoch):
        #t3 = time.time()        
        
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
        
        rootsListUTC = self.seconds_2_utc(epoch, rootsList)
        
        for i in range(0, len(rootsListUTC)):
            rootsListUTC[i] = self.datetime_2_timestamp(rootsListUTC[i]) 
       
        #print(len(rootsListUTC))        
        
        if(VF[0] > 0):  #then we are already in visibility range, so t0 to the first root is an interval
            #timingWindow.append([0, rootsList[0]])
            timingWindow.append([epoch, rootsListUTC[0]])
            startIndex = 1
        elif(VF[0] < 0):  #then we are not in visibility range, so the first root will be the start of an interval
            #timingWindow.append([rootsList[0], rootsList[1]])
            timingWindow.append([rootsListUTC[0], rootsListUTC[1]])
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
                #timingWindow = [0, rootsList[endrootIndex]]
                timingWindow = [epoch, rootsListUTC[endrootIndex]]
                startIndex = endrootIndex+1
            else:
                #timingWindow = [rootsList[endrootIndex], rootsList[endrootIndex+1]]
                timingWindow = [rootsListUTC[endrootIndex], rootsListUTC[endrootIndex+1]]
                startIndex = endrootIndex+2
        
#        endIndex = len(rootsList)+1
        endIndex = len(rootsListUTC) - 1
       
        for index in range(startIndex, endIndex, 2):
            timingWindow.append([rootsListUTC[index], rootsListUTC[index+1]])

        
        #Check to see if the roots form a closed set.
        
        if(range(startIndex, endIndex)[-1] == (endIndex-2)):          
            if(rootsList[-1] != times[-1]):
                timingWindow.append([rootsListUTC[-1], self.datetime_2_timestamp(self.one_time_to_utc(epoch, times[-1]))])
                #timingWindow.append([rootsListUTC[-1], times[-1]])
        
        #if(self.verbose):
        #for window in timingWindow:
         #   print(window)
        
        #print("Time for creating timing windows = "+str(time.time() - t3))        
        
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
        
        
    def satellite_viewing_cone(self, dataECEF, posECEF, posGEO):
        r_sat = dataECEF[:,0:3]  
        r_site = posECEF
        v_sat = dataECEF[:,3:6]   
        lat = np.radians(posGEO[0])
        lon = posGEO[1]
        theta = 0 #Not used yet.
        m = 0
        
        m_r_site = np.sqrt(np.sum(r_site*r_site, axis=1))      
        m_r_sat = np.sqrt(np.sum(r_sat*r_sat, axis=1))    
        m_q = np.max(m_r_sat)        
        
        gamma = theta + np.arcsin( (m_r_site*np.sin((np.pi/2.0)+theta))/(m_q) )

        p = np.cross(r_sat, v_sat)
        p_m = np.sqrt(np.sum(p*p, axis=1))
        p_unit = p/p_m[:,None]
        
        tInOut = (np.arcsin( (np.cos(gamma) - p_unit[:,2]*np.sin(lat))/
            (np.sqrt(np.power(p_unit[:,0],2.0) + np.power(p_unit[:,1],2.0))*np.cos(lat)) )
                - lon - np.arctan(p_unit[:,0]/p_unit[:,1] + 2*np.pi*m))
        
        return np.array(tInOut)
        
    def seconds_2_utc(self, epoch, times):
        utc_times=[]
        
        for i in range(0,len(times)):
            deltat = datetime.timedelta(seconds=times[i])
            utc_times.append(epoch + deltat)
        
        return utc_times
        
        
    def one_time_to_utc(self, epoch, time):
        deltat = datetime.timedelta(seconds=time)
        
        return epoch + deltat
        
        
    def calc_epoch(self, jdate):
        days = jdate - self.RataDieJDtime - 1
        setup_delta = datetime.timedelta(days=days)
        epoch = datetime.datetime(1,1,1) + setup_delta
        
        return epoch
        
    
    def check_time_intersection(self, m_start, m_end, d_start, d_end, times):
        if(m_start < d_start):
            start = d_start
            start_index = 0     
        else:
            start = m_start
            start_jump = (start - d_start).total_seconds()
            start_index = self.binary_List_Search(times, start_jump)
        if(m_end > d_end):
            end = d_end
            end_index = -1
        else:
            end = m_end
            end_jump = (d_end - end).total_seconds()
            end_index = self.binary_List_Search(times, end_jump)

        if(self.verbose):
            print(start)
            print(end)
        
        return [start_index, end_index]
        
        
    def datetime_2_timestamp(self, dt):
        year = str(dt.date().year)
        month = self.resize_2(str(dt.date().month))
        day = self.resize_2(str(dt.date().day))
        hour = self.resize_2(str(dt.time().hour))
        minute = self.resize_2(str(dt.time().minute))
        second = self.resize_2(str(dt.time().second))
        milsec = self.resize_3(str(dt.time().microsecond/1000.0))        
        
        dateString = year+month+day
        timeString = hour+":"+minute+":"+second+"."+milsec
        
        datetimeString = dateString+"T"+timeString
        
        return datetimeString
        
    
    def resize_2(self, myString):
        if(len(myString) == 1):
            myString = "0"+myString
        elif(len(myString) > 2):
            myString = myString[0:2]
        
        return myString
    
    
    def resize_3(self, myString):
        if(len(myString) == 1):
            myString = "00"+myString
        elif(len(myString) == 2):
            myString = "0"+myString
        elif(len(myString) > 3):
            myString = myString[0:3]
        
        return myString