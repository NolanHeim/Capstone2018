# -*- coding: utf-8 -*-
#
# Sensor.py
#
# Sensor Calculations
#   Includes the algorithms to determine the sensor viewing projection on
#   the surface of the Earth. Also includes the methods to determine the
#   intersection between two polygons (i.e. the AOI and viewing area)
#
# Initial Creation Date: 11/17/2018
#
# Written by Jordan Jones and Nolan Heim
#

import numpy as np
import math
from Calculator import *
import time as tm
from sympy.geometry import *

class Sensor:    
    def __init__(self):
        self.equitorialRadius = 6378137.0 #m
        self.polarRadius = 6356752.3 #m
        self.EARTH_RADIUS = 6371000.0 #m

    #Determines the intersection between the Area of interest (AOI) and 
    # satellite viewing area
    def sensors_intersection(self, satelliteData, sensors, mission, Npoly, delta_t):
        sensorType = mission.get_sensor_type()
        time = np.array(satelliteData[:,0])
        booleanTimes = False
        #Determine the sensor type
        if(sensorType == ''):
            #Cycle through all sensors
            for sensor in sensors:
                booleanTimeSensor = self.single_sensor_intersection(satelliteData, sensor, Npoly, delta_t)
                booleanTimes = booleanTimes | booleanTimeSensor
        else:
            for sensor in sensors:
                if(sensor['Imaging Type'] == sensorType):
                    booleanTimeSensor = single_sensor_intersection(satelliteData, sensor, Npoly, delta_t)
                    booleanTimes = booleanTimes | booleanTimeSensor
        
        booleanTimes = np.array(booleanTimes)
        listOfTimeWindows = time[booleanTimes]     

        boundCondition = np.diff(listOfTimeWindows) > delta_t
        rightBounds = (np.argwhere(boundCondition).T)[0]
        if(len(rightBounds) == 0): #no right bounds
            return np.array([])
        
        if(rightBounds[-1] != len(listOfTimeWindows)-1):
            rightBounds = np.append(rightBounds, len(listOfTimeWindows)-1)
        if(rightBounds[0] == 0):
            rightBounds = rightBounds[1:]
            
        leftBounds = np.array(rightBounds[:-1]) + 1
        leftBounds = np.insert(leftBounds, 0, 0)
            
        #Using the left & right bounds for each point, construct the timing windows. 
        timeWindows = []
        for bound in range(0, len(rightBounds)):
            timeWindows.append([listOfTimeWindows[leftBounds[bound]], 
                                listOfTimeWindows[rightBounds[bound]]])
 
        timeWindows = np.array(timeWindows)
        
        return timeWindows

    #Assumed Npoly of the form: [point1, point2, point3, point4] in [lat,long]
    #Determines the intersection between the satellite and
    #the N-D polygon
    #Note: Due to the sympy intersection library this method is not optimized.
    def single_sensor_intersection(self, satelliteData, sensor, Npoly, delta_t):                
        viewingType = sensor['Angle Dependent']        
        time = satelliteData[:,0]
        
        #Determine if sensor has independent or dependent range of motion angles
        if(viewingType == "False"):
            #A rectangular viewing area
            Npoly = [Point(i) for i in Npoly]        
            AOI = Polygon(*Npoly)
            [rectangleECEF, rArea] = self.viewing_rectangle(satelliteData, sensor)
            timeIndex = 0
            timeWindowIndex = []           
            for t in time:
                viewingRectangle = np.array([
                self.ecef_to_geo_on_surface(rectangleECEF[0][timeIndex]),
                self.ecef_to_geo_on_surface(rectangleECEF[1][timeIndex]),
                self.ecef_to_geo_on_surface(rectangleECEF[2][timeIndex]),
                self.ecef_to_geo_on_surface(rectangleECEF[3][timeIndex])])
                
                pointViewingArea = []
                for point in viewingRectangle:
                    pointViewingArea.append(Point(point[0][0],point[0][1]))
            
                
                VA = Polygon(*pointViewingArea)                         
                if(intersection(VA,AOI) != []):
                    #Intersection between AOI and Viewing Area
                    timeWindowIndex.append(timeIndex)
                timeIndex = timeIndex + 1
                
            timeWindows = time[timeWindowIndex]
            #Take the intersection between the total time and the subset
            #Results in a boolean array
            booleanTimes = np.in1d(time,timeWindows)
        else:
            booleanTimes = False                
            Npoly = self.geo_to_ecef(Npoly[0],Npoly[1],0.0*Npoly[0])
            #A circular viewing area (dependent range of motion angle case)
            [centroid, radius] = self.viewing_circle(satelliteData, sensor)
            for pointn in Npoly:
                d_rn = centroid - pointn
                m_d_rn = np.sqrt(np.sum(d_rn*d_rn, axis=1))
                radiusFunction = m_d_rn - radius
                booleanTimes = booleanTimes | (radiusFunction < 0)
        
        return booleanTimes

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

    #Determines the conversion from ECEF to geo assuming an
    #altitude of 0 (on Earth's surface)
    def ecef_to_geo_on_surface(self, position):
        a = self.equitorialRadius
        b = self.polarRadius

        x = position[0]
        y = position[1]
        z = position[2]        
        
        lam = np.arctan2(y,x)

        e = np.sqrt(np.power(a,2.0)-np.power(b,2.0))/a
        p = np.sqrt(np.power(x,2.0)+np.power(y,2.0))
        alt = 0.0
        phi = math.atan2(z,p*(1.0-np.power(e,2.0)))

        lon = np.degrees(lam)
        lat = np.degrees(phi)

        geo = np.column_stack([lat, lon, alt])
        
        return geo
     
    #based off the second response here https://gis.stackexchange.com/questions/265909/converting-from-ecef-to-geodetic-coordinates
    def ecef_to_geo(self, position):
        a = self.equitorialRadius
        b = self.polarRadius

        x = position[0]
        y = position[1]
        z = position[2]        
        
        lam = np.arctan2(y,x)

        e = np.sqrt(np.power(a,2.0)-np.power(b,2.0))/a
        p = np.sqrt(np.power(x,2.0)+np.power(y,2.0))
        
        # have to iterate to find alt, with an initial guess of 0
        alt_0 = 0.0
        phi = math.atan2(z,p*(1.0-np.power(e,2.0)))
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        N = np.power(a,2.0)/np.sqrt(np.power(a*cos_phi,2.0)+np.power(b*sin_phi,2.0))
        alt = p/cos_phi - N
        while abs(alt-alt_0) > 1.0e-3:
            alt_0 = alt
            phi = np.arctan2(z,p*(1.0-np.power(e,2.0)*N/(N+alt)))
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
        
    
    #Returns the for corners of the viewing rectangle in ECEF
    def viewing_rectangle(self, satelliteData, sensor):
        dphi_h = sensor["Rotation"][0]*(np.pi)/180
        dphi_v = sensor["Rotation"][1]*(np.pi)/180        
        
        [centroid, d_centroid, h_unit, v_unit] = self.determineCentroid(satelliteData, sensor)

        #Need to construct the viewing rectangle
        #Determine the projection of dphi onto plane
        deltah = d_centroid*np.tan(dphi_h)
        deltav = d_centroid*np.tan(dphi_v)
        
        point1 = centroid - h_unit*deltah - v_unit*deltav
        point2 = centroid - h_unit*deltah + v_unit*deltav
        point3 = centroid + h_unit*deltah + v_unit*deltav
        point4 = centroid + h_unit*deltah - v_unit*deltav
        
        ## Rect. Corner Points orientation
        # |2 3|
        # |1 4|
        ##
        rectangle = [point1, point2, point3, point4]
        rArea = (2*deltah)*(2*deltav)        
        
        return [rectangle, rArea]
    
    #Determines the centroid point
    #This represents the point on Earth's surface where the central viewing
    #axis of the satellite intersects.
    def determineCentroid(self, satelliteData, sensor):     
        phi_h = sensor["Orientation"][0]*(np.pi)/180
        phi_v = sensor["Orientation"][1]*(np.pi)/180
        
        r_sat = satelliteData[:,1:4]
        #Assume a constant distance from the center of the earth
        m_r_sat = np.sqrt(np.power(r_sat[:,0], 2.0) + np.power(r_sat[:,1], 2.0) + np.power(r_sat[:,2], 2.0))
        m_r_sat_avg = np.average(m_r_sat)
        v_sat = satelliteData[:,4:7]
        
        a = m_r_sat_avg - self.EARTH_RADIUS 
        #Part 1: Determine the distance along the centroid vector to the earth's surface
        #Step 1: Determine Alpha
        dv = a*np.tan(phi_v)
        dh = a*np.tan(phi_h)
        alpha = np.sqrt(np.power(dv, 2.0) + np.power(dh, 2.0))/a
        #Step 2: Determine the total length of the centroid distance (l)
        theta = np.arcsin(self.EARTH_RADIUS/m_r_sat_avg)
        lPrime = m_r_sat_avg*np.cos(theta)
        l = lPrime/np.cos(theta - alpha)
        #Step 3: Determine distance along centroid vector (dCent)
        k = np.sqrt(np.power(l, 2.0) + np.power(m_r_sat_avg, 2.0) - 2.0*l*m_r_sat_avg*np.cos(alpha))
        beta = np.pi - ((np.pi/2) - (theta - alpha))
        gamma = np.arcsin((k*np.sin(beta))/(self.EARTH_RADIUS))
        delta = np.pi - beta - gamma
        c = np.sqrt(np.power(k, 2.0) + np.power(self.EARTH_RADIUS, 2.0) - 2*k*self.EARTH_RADIUS*np.cos(delta))
        dCent = (l-c)
        
        #Determine the orthogonal unit vector (h_unit)
        m_r_sat = np.sqrt(np.sum(r_sat*r_sat, axis=1))
        r_unit_sat = r_sat/m_r_sat[:,None]
        m_v_sat = np.sqrt(np.sum(v_sat*v_sat, axis=1))
        v_unit = v_sat/m_v_sat[:,None]
        
        h_vec = np.cross(r_unit_sat, v_unit)
        m_h_vec = np.sqrt(np.sum(h_vec*h_vec, axis=1))
        h_unit = h_vec/m_h_vec[:,None]     
        
        #Part 2: Determine the unit vector for the centroid point from satellite
        rotVCent = self.computeRotation(phi_v, phi_h, h_unit, v_unit, (-1)*r_sat)
        m_rotVCent = np.sqrt(np.sum(rotVCent*rotVCent, axis=1))
        unit_cent_sat = rotVCent/m_rotVCent[:,None]
        
        #Part 3: Determine the vector for the centroid point from the origin of ECEF
        cent_sat = unit_cent_sat*dCent
        #Determine the point on the earth corresponding to the centroid in ECEF
        cent_ecef = r_sat + cent_sat
        
        return [cent_ecef, dCent, h_unit, v_unit]
    
    #Computes the rotation for the centroid calculation
    def computeRotation(self, phi_v, phi_h, h_unit, v_unit, vec):
            
        vRot = self.arbitraryRotationMatrix(phi_v, h_unit)
        hRot = self.arbitraryRotationMatrix(phi_h, v_unit)
        
        #Define the new rotational Matrix (vRot*hRot)
        fRot = {}
        #Row 0         
        fRot[0] = vRot[0]*hRot[0] + vRot[1]*hRot[3] + vRot[2]*hRot[6]
        fRot[1] = vRot[0]*hRot[1] + vRot[1]*hRot[4] + vRot[2]*hRot[7]
        fRot[2] = vRot[0]*hRot[2] + vRot[1]*hRot[5] + vRot[2]*hRot[8]
        #Row 1
        fRot[3] = vRot[3]*hRot[0] + vRot[4]*hRot[3] + vRot[5]*hRot[6]
        fRot[4] = vRot[3]*hRot[1] + vRot[4]*hRot[4] + vRot[5]*hRot[7]
        fRot[5] = vRot[3]*hRot[2] + vRot[4]*hRot[5] + vRot[5]*hRot[8]
        #Row 2
        fRot[6] = vRot[6]*hRot[0] + vRot[7]*hRot[3] + vRot[8]*hRot[6]
        fRot[7] = vRot[6]*hRot[1] + vRot[7]*hRot[4] + vRot[8]*hRot[7]
        fRot[8] = vRot[6]*hRot[2] + vRot[7]*hRot[5] + vRot[8]*hRot[8]
        
        #Comput the unit vector for the given vector (vec)
        rotVec = {}
        rotVec[0] = vec[:,0]*fRot[0] + vec[:,1]*fRot[1] + vec[:,2]*fRot[2]
        rotVec[1] = vec[:,0]*fRot[3] + vec[:,1]*fRot[4] + vec[:,2]*fRot[5]
        rotVec[2] = vec[:,0]*fRot[6] + vec[:,1]*fRot[7] + vec[:,2]*fRot[8]
        
        rotVecMatrix = np.column_stack((rotVec[0], rotVec[1], rotVec[2]))

        return rotVecMatrix
            
    #Returns the elements for a rotational matrix for
    # a rotational about an arbitrary axis (u) and angle (a)
    #   | 0 1 2 |
    #   | 3 4 5 |
    #   | 6 7 8 |
    def arbitraryRotationMatrix(self, a, u):
        ux = u[:,0]
        uy = u[:,1]
        uz = u[:,2]
        rotationMatrix = {}
        #Row 0
        rotationMatrix[0] = np.cos(a) + np.power(ux, 2.0)*(1-np.cos(a))
        rotationMatrix[1] = ux*uy*(1-np.cos(a)) - uz*np.sin(a)
        rotationMatrix[2] = ux*uz*(1-np.cos(a)) + uy*np.sin(a)
        #Row 1
        rotationMatrix[3] = uy*ux*(1-np.cos(a)) + uz*np.sin(a)
        rotationMatrix[4] = np.cos(a) + np.power(uy, 2.0)*(1-np.cos(a))
        rotationMatrix[5] = uy*uz*(1-np.cos(a)) - ux*np.sin(a)
        #Row 2       
        rotationMatrix[6] = uz*ux*(1-np.cos(a)) - uy*np.sin(a)
        rotationMatrix[7] = uz*uy*(1-np.cos(a)) + ux*np.sin(a)
        rotationMatrix[8] = np.cos(a) + np.power(uz, 2.0)*(1-np.cos(a))
                
        return rotationMatrix
        
    def get_psi(self, phi, altitude):
        theta = np.asin((altitude + self.r_earth)*(np.sin(phi)/self.r_earth))
        return (np.pi - theta - phi)
        
    def viewing_circle(self, satelliteData, sensor):
        #Get the centroid of the viewing circle on earth
        dphi_c = sensor['Rotation']*(np.pi)/180
        [centroid, dCentroid, h_unit, v_unit] = self.determineCentroid(satelliteData, sensor)
        
        #Now determine the radius of the circle
        radius = dCentroid*np.tan(dphi_c)

        return [centroid, radius]        
        
        