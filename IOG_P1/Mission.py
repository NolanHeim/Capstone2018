# -*- coding: utf-8 -*-
#
# Mission.py
#
# (PHASE 1)
#
# Represents a mission, contains all needed parameters for calculator. 
#
# Initial Creation Date: 09/28/2018
#
# Written by Jordan Jones and Nolan Heim
#
import datetime

class Mission:
    
    #set fields based on parameters
    def __init__(self, targetCoordinates, name, 
                 intervalStart, intervalEnd, idsToConsider):
        self.targetCoordinates = targetCoordinates
        self.name = name
        self.intervalStart = intervalStart
        self.intervalEnd = intervalEnd
        self.idsToConsider = idsToConsider

        #Parse the start and end times here.
        #Save as datetime objects.
        self.parse_date_time()
        

    #Handles parsing the intervalStart/intervalEnd strings as specified in the 
    #REST API spec. Numbers are hard coded based on that string, as we assume
    #the format willl not change.
    def parse_date_time(self):
        startYear = int(self.intervalStart[0:4])
        startMonth = int(self.intervalStart[4:6])
        startDay = int(self.intervalStart[6:8])
        startHour = int(self.intervalStart[9:11])
        startMinute = int(self.intervalStart[12:14])
        startSec = int(self.intervalStart[15:17])
        
        self.startDateTime = datetime.datetime(startYear, startMonth, startDay,
                                               startHour, startMinute, startSec)
        
        endYear = int(self.intervalEnd[0:4])
        endMonth = int(self.intervalEnd[4:6])
        endDay = int(self.intervalEnd[6:8])
        endHour = int(self.intervalEnd[9:11])
        endMinute = int(self.intervalEnd[12:14])
        endSec = int(self.intervalEnd[15:17])
        

        self.endDateTime = datetime.datetime(endYear, endMonth, endDay,
                                             endHour, endMinute, endSec)    
        
        
    def get_name(self):
        return self.name
    
    
    def get_coordinates(self):
        self.targetCoordinates.append(0.0) #Accounting for altitude
        return self.targetCoordinates
        
        
    def get_interval_start_time(self):
        return self.startDateTime
        
    
    def get_interval_end_time(self):
        return self.endDateTime

        
    def get_ids_to_consider(self):
        return self.idsToConsider
        
        
    #Debugging function only    
    def display_parameters(self):
        print(self.get_coordinates())
        print(self.get_interval_start_time())
        print("Mission ID "+str(self.get_name()))
        print(self.get_ids_to_consider())
        
        
    #Used to ensure all required parameters have values    
    def check_params(self):
        if (self.targetCoordinates == [] or self.name == "" or self.intervalEnd == ""
        or self.intervalStart == ""):
            return False
        else:
            return True