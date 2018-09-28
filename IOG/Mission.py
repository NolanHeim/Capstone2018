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

class Mission:
    
    def __init__(self, targetCoordinates, name, sensorType, illumDir, illumThresh, intervalStart, intervalEnd):
        self.targetCoordinates = targetCoordinates
        self.name = name
        self.sensorType = sensorType
        self.illumDir = illumDir
        self.illumThresh = illumThresh
        self.intervalStart = intervalStart
        self.intervalEnd = intervalEnd
    
    
    def get_name(self):
        return self.name
    
    
    def get_coordinates(self):
        return self.targetCoordinates
        
    
    def get_sensor_type(self):
        return self.sensorType
        
    
    def get_illumination_direction(self):
        return self.illumDir
        
        
    def get_illumination_threshold(self):
        return self.illumThresh
        
        
    def get_interval_start_time(self):
        return self.intervalStart
        
    
    def get_interval_end_time(self):
        return self.intervalEnd