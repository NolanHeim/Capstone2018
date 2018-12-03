# -*- coding: utf-8 -*-
#
# Satellite.py
#
# Represents the data for a single satellite
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#

import numpy as np

class Satellite:  
    def __init__(self, extraInfo):
        self.epoch_index = 0       
        self.uuid_index = 1
        self.filename_index = 2
        self.parsed_filename_index = 3
        
        self.extraInfo = self.sort_extra_info(extraInfo)
        
        self.sensors = []
        self.dataMatrix = 0

        
    def sort_extra_info(self, extraInfo):
        sortedInfo = {}
        sortedInfo['epoch'] = extraInfo[self.epoch_index]
        sortedInfo['uuid'] = extraInfo[self.uuid_index]
        sortedInfo['filename'] = extraInfo[self.filename_index]
        sortedInfo['parsed_filename'] = extraInfo[self.parsed_filename_index]
        return sortedInfo

    
    def set_sensors(self, sensor_model): #sensor_model is a dictionary containing all sensor
        for sensor in sensor_model:
            self.sensors.append(sensor_model[sensor])
            
            
    def set_data_matrix(self, dataMatrix):
        self.dataMatrix = dataMatrix
        print("we have been set")
        print(np.array(dataMatrix).shape)
    
    
    def get_data_matrix(self):
        return np.array(self.dataMatrix)

        
    def get_epoch(self):
        return self.extraInfo['epoch']

    
    def get_uuid(self):
        return self.extraInfo['uuid']

        
    def get_satellite_name(self):
        return self.extraInfo['filename']


    def get_parsed_orbit_name(self):
        return self.extraInfo['parsed_filename']
        
        
    def get_sensors(self):
        return self.sensors