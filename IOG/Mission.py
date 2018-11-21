# -*- coding: utf-8 -*-
#
# Mission.py
#
# Represents a mission, contains all needed parameters for calculator. 
#
# Initial Creation Date: 09/28/2018
#
# Written by Jordan Jones and Nolan Heim
#
import numpy as np
import datetime

class Mission:
    
    def __init__(self, targetCoordinates, name, sensorType, illumDir, illumThresh, 
                 intervalStart, intervalEnd, idsToConsider):
        self.targetCoordinates = targetCoordinates
        self.name = name
        self.sensorType = sensorType
        self.illumDir = illumDir
        self.illumThresh = illumThresh
        self.intervalStart = intervalStart
        self.intervalEnd = intervalEnd
        self.idsToConsider = idsToConsider

        #Parse the start and end times here.
        #Save as datetimes.
        self.parse_date_time()
        

    def parse_date_time(self):
        startYear = int(self.intervalStart[0:4])
        startMonth = int(self.intervalStart[4:6])
        startDay = int(self.intervalStart[6:8])
        startHour = int(self.intervalStart[9:11])
        startMinute = int(self.intervalStart[12:14])
        startSec = int(self.intervalStart[15:17])
        #startMilSec = int(self.intervalStart[18:-1])        
        
        self.startDateTime = datetime.datetime(startYear, startMonth, startDay,
                                               startHour, startMinute, startSec)
        
        endYear = int(self.intervalStart[0:3])
        endMonth = int(self.intervalStart[4:5])
        endDay = int(self.intervalStart[6:7])
        endHour = int(self.intervalStart[9:10])
        endMinute = int(self.intervalStart[12:13])
        endSec = int(self.intervalStart[15:16])
        #endMilSec = int(self.intervalStart[18:-1])        
        
        
        self.endDateTime = datetime.datetime(endYear, endMonth, endDay,
                                             endHour, endMinute, endSec)
        
        
        #startTime = self.intervalStart.split()
        #startTime[1] = startTime[1].replace("GMT", "")
        #endTime = self.intervalEnd.split()
        #endTime[1] = endTime[1].replace("GMT","")
        #start_mmddyyyy = startTime[0].split('.')
        #start_hhmmss = startTime[1].split(':')
        #end_mmddyyyy = endTime[0].split('.')
        #end_hhmmss = endTime[1].split(':')
        
        #start_mmddyyyy = [int(i) for i in start_mmddyyyy]        
        #start_hhmmss = [int(i) for i in start_hhmmss]
        #end_mmddyyyy = [int(i) for i in end_mmddyyyy]        
        #end_hhmmss = [int(i) for i in end_hhmmss]        

        #self.startDateTime = datetime.datetime(start_mmddyyyy[2], start_mmddyyyy[0], start_mmddyyyy[1], 
        #                                       start_hhmmss[0], start_hhmmss[1], start_hhmmss[2])
        #self.endDateTime = datetime.datetime(end_mmddyyyy[2], end_mmddyyyy[0], end_mmddyyyy[1], 
        #                                       end_hhmmss[0], end_hhmmss[1], end_hhmmss[2])
        

    def get_name(self):
        return self.name
    
    
    def get_coordinates(self):
        self.targetCoordinates.append(0.0) #Accounting for altitude
        return self.targetCoordinates
        
    
    def get_sensor_type(self):
        return self.sensorType
        
    
    def get_illumination_direction(self):
        return self.illumDir
        
        
    def get_illumination_threshold(self):
        return self.illumThresh
        
        
    def get_interval_start_time(self):
        return self.startDateTime
        
    
    def get_interval_end_time(self):
        return self.endDateTime
        
    def get_ids_to_consider(self):
        return self.idsToConsider
        
    def display_parameters(self):
        print(self.get_coordinates())
        print(self.get_interval_start_time())
        print("Mission ID "+str(self.get_name()))
        print(self.get_ids_to_consider())
        
    def check_params(self):
        if (self.targetCoordinates == [] or self.name == "" or self.intervalEnd == ""
        or self.intervalStart == "" or self.sensorType == ""):
            return False
        else:
            return True