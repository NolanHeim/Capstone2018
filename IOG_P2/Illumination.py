# -*- coding: utf-8 -*-
#
# Illumination.py
#
# Represents the data for a single satellite
#
# Initial Creation Date: 10/27/2018
#
# Written by Jordan Jones and Nolan Heim
#
import numpy as np

class Illumination:
    #TODO insert class functions and parameters    
    def __init__(self):
        self.start = True
    
    #time is a datetime object.
    #Lat/Lon in degrees
    #https://support.microsoft.com/en-ca/help/214019/method-to-determine-whether-a-year-is-a-leap-year
    def computeSolarAngles(self, Lat, Lon, time):
        #Determine if year is a leap year
        year = time.year
        if(year % 4 == 0):
            if(year % 100 == 0):
                if(year % 400 == 0):
                    daysInYear = 366.0;
                else:                        
                    daysInYear = 365.0;
            else:
                daysInYear = 366.0;
        else:
            daysInYear = 365.0;
        
        #Compute fractional year (y)
        y = ((2*np.pi)/(daysInYear))*(time.timetuple().tm_yday - 1.0 + (time.hour)/(24.0))
        eqtime = 229.18*(0.000075 + 0.001868*np.cos(y) - 0.032077*np.sin(y) - 
                    0.014615*np.cos(2.0*y) - 0.040849*np.sin(2.0*y) )
        decl = (0.006918 - 0.399912*np.cos(y) + 0.070257*np.sin(y) - 
                    0.006758*np.cos(2.0*y) + 0.000907*np.sin(2.0*y) - 
                    0.002697*np.cos(3.0*y) + 0.00148*np.sin(3.0*y))
        time_offset = eqtime + 4.0*Lon
        tst = time.hour*60.0 + time.minute + time.second/60.0 + time_offset
        ha = (tst/4.0) - 180.0
        #print(ha)
        radLat = np.radians(Lat)
        #print(radLat)
        solarZenith = np.arccos(np.sin(radLat)*np.sin(decl) +
                        np.cos(radLat)*np.cos(decl)*np.cos(np.radians(ha)))
        solarElevation = (np.pi/2.0) - solarZenith
        if(time.hour > 12):
           solarAzimuth = 2.0*np.pi - np.arccos((-1.0)*(np.sin(radLat)*np.cos(solarZenith) - np.sin(decl))/
                            (np.cos(radLat)*np.sin(solarZenith)))
        else:
            solarAzimuth = np.arccos((-1.0)*(np.sin(radLat)*np.cos(solarZenith) - np.sin(decl))/
                            (np.cos(radLat)*np.sin(solarZenith)))
        
        #Change to degrees
        solarInclination = 90 + solarElevation*(180.0/np.pi)

        return solarInclination     
        
    #Takes solarElevation as a numpy array as input.
    def computeSolarIrradiation(self, solarElevation):
        solarZenith = (np.pi/2.0) - solarElevation
        
        #
        irradiation = 910*np.cos(solarZenith) - 30;
        
        negativeFlag = irradiation < 0.0
        irradiation[negativeFlag] = 0
        #W/(m^2)
        return irradiation
