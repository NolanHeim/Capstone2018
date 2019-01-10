# -*- coding: utf-8 -*-
#
# Satellite.py
#
# (PHASE 2)
#
# Represents the data for a single satellite in a Constellation.
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#

import numpy as np

class Satellite:  
    def __init__(self, extraInfo):
        #represent where in the extraInfo each piece of information can be found (see Parser for order)
        self.epoch_index = 0       
        self.uuid_index = 1
        self.filename_index = 2
        self.parsed_filename_index = 3
        
        self.extraInfo = self.sort_extra_info(extraInfo)
        
        self.sensors = []
        self.dataMatrix = 0


    # Converts the extra info from a list to a dictionary
    def sort_extra_info(self, extraInfo):
        sortedInfo = {}
        sortedInfo['epoch'] = extraInfo[self.epoch_index]
        sortedInfo['uuid'] = extraInfo[self.uuid_index]
        sortedInfo['filename'] = extraInfo[self.filename_index]
        sortedInfo['parsed_filename'] = extraInfo[self.parsed_filename_index]
        return sortedInfo

    
    def set_sensors(self, sensor_model): #sensor_model is a dictionary containing all sensor objects
        for sensor in sensor_model:
            self.sensors.append(sensor_model[sensor])
            
    
    def set_data_matrix(self, dataMatrix):
        self.dataMatrix = dataMatrix
        
    
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