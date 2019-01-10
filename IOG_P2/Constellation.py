# -*- coding: utf-8 -*-
#
# Constellation.py
#
# (PHASE 2)
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
        
    
    # Loads saved satellite data into each satellite object. uuids is the list of
    # uuids to consider, specified by the Mission.
    def add_satellite_data(self, uuids):
        if(uuids == []): #then load data for ALL available satellites
            for uuid in self.constellation:
                satellite = self.constellation[uuid]
                dataMatrix = self.load_data_matrix(satellite.get_parsed_orbit_name())
                sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                satellite.set_data_matrix(dataMatrix)
                satellite.set_sensors(sensorModel)

        else: #only load data for specified satellites 
            for uuid in uuids:
                if uuid in self.constellation:
                    satellite = self.constellation[uuid]
                    dataMatrix = self.load_data_matrix(satellite.get_parsed_orbit_name())
                    sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                    satellite.set_data_matrix(dataMatrix)
                    satellite.set_sensors(sensorModel)
                else:
                    return 'ERROR: UUID is not valid'


    # Returns the Satellite object corresponding to the given UUID
    def get_satellite(self, UUID):
        return self.constellation[UUID]

        
    # Gets a list containing all the Satellite instantiations in this Constellation
    def get_satellite_list(self):
        satList = [];
        for uuid in self.constellation:
            satList.append(self.constellation[uuid])
        return satList


    # Instantiates Satellite objects for each extraInfo file
    def construct_satellites(self):
        for ei in self.extraInfoMatrix:
            row = [ei["jDate"], ei["satelliteID"], ei["satelliteName"], ei["parsed_dataFile"]]
            sat = Satellite(row)
            self.constellation[sat.get_uuid] = sat

        
   # Creates a memmap to read parsed satellite data. Note that this function assumes the parser
   # has already been ran by itself.
    def load_data_matrix(self, filename):
        file_with_extension = filename
        if file_with_extension in os.listdir(self.parsed_datapath):
            matrix = np.load(self.parsed_datapath+file_with_extension, mmap_mode='r')
        else:
            matrix = 0
            
        return np.array(matrix)


    # Loads sensors fo a given satellite
    def load_sensor_model(self, satName):
        filename = 'sensor_parsed_' + satName + '.json'
        with open(self.parsed_datapath + filename, "r") as sensor_json:
                sensorModel = json.load(sensor_json)       
        return sensorModel


    # Obtains all the extraInfo all at once and concatenates it all into a list of dictionaries
    def load_extra_info(self):
        ei_Matrices = []

        for filename in os.listdir(self.parsed_datapath):
            if(filename[0] == 'e'): #begins with e => must be extra info (see spec.)           
                with open(self.parsed_datapath + filename, "r") as extra_info_json:
                    extraInfoDict = json.load(extra_info_json)       
        
                ei_Matrices.append(extraInfoDict)
            
        return np.array(ei_Matrices)    