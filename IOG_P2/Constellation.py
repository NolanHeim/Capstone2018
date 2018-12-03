# -*- coding: utf-8 -*-
#
# Constellation.py
#
# Represents the group of all satellites for which data is provided.
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
from Satellite import *
import numpy as np
import os
import json

class Constellation:    
    def __init__(self, parsed_datapath):
        self.constellation = {};
        self.parsed_datapath = parsed_datapath
        self.extraInfoMatrix = self.load_extra_info()
        self.construct_satellites()
        
    
    def add_satellite_data(self, uuids):
        #then load data for all available satellites
        if(uuids == []): 
            for uuid in self.constellation:
                satellite = self.constellation[uuid]
                dataMatrix = self.load_data_matrix(satellite.get_parsed_orbit_name())
                sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                satellite.set_data_matrix(dataMatrix)
                satellite.set_sensors(sensorModel)

        #only load data for specified satellites
        else: 
            for uuid in uuids:
                if uuid in self.constellation:
                    satellite = self.constellation[uuid]
                    dataMatrix = self.load_data_matrix(satellite.get_parsed_orbit_name())
                    sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                    satellite.set_data_matrix(dataMatrix)
                    satellite.set_sensors(sensorModel)
                else:
                    return 'ERROR: UUID is not valid'

            
    def get_satellite(self, UUID):
        return self.constellation[UUID]
        
        
    def get_satellite_list(self):
        satList = [];
        for uuid in self.constellation:
            satList.append(self.constellation[uuid])
        return satList

    
    def construct_satellites(self):
        for ei in self.extraInfoMatrix:
            row = [ei["jDate"], ei["satelliteID"], ei["satelliteName"], ei["parsed_dataFile"]]
            sat = Satellite(row)
            self.constellation[sat.get_uuid] = sat

        
   # Creates a memmap to read parsed satellite data. Note that this function assumes the parser
   # has already been ran by itself.
    def load_data_matrix(self, filename):
        file_with_extension = filename
        print(file_with_extension)
        if file_with_extension in os.listdir(self.parsed_datapath):
            matrix = np.load(self.parsed_datapath+file_with_extension, mmap_mode='r')
            print("we found the orbit file")
        else:
            matrix = 0
            
        return np.array(matrix)

        
    def load_sensor_model(self, satName):
        filename = 'sensor_parsed_' + satName + '.json'
        with open(self.parsed_datapath + filename, "r") as sensor_json:
                sensorModel = json.load(sensor_json)       
        return sensorModel

        
    def load_extra_info(self):
        ei_Matrices = []

        for filename in os.listdir(self.parsed_datapath):
            if(filename[0] == 'e'):            
                with open(self.parsed_datapath + filename, "r") as extra_info_json:
                    extraInfoDict = json.load(extra_info_json)       
        
                ei_Matrices.append(extraInfoDict)
            
        return np.array(ei_Matrices)    