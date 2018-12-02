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
                dataMatrix = self.load_data_matrix(satellite.get_satellite_name())
                sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                satellite.set_data_matrix(dataMatrix)
                satellite.set_sensors(sensorModel)

        #only load data for specified satellites
        else: 
            for uuid in uuids:
                if uuid in self.constellation:
                    satellite = self.constellation[uuid]
                    dataMatrix = self.load_data_matrix(satellite.get_satellite_name())
                    sensorModel = self.load_sensor_model(satellite.get_satellite_name())
                    satellite.set_data_matrix(dataMatrix)
                    satellite.set_sensors(sensorModel)
                else:
                    return 'ERROR: UUID is not valid'

            
    def get_satellite(self, UUID):
        return self.constellation[UUID]
        
        
    def get_satellite_list(self):
        satList = [];
        for satellite in self.constellation:
            satList.append(satellite)
        return satList

    
    def construct_satellites(self):
        for row in self.extraInfoMatrix:
            sat = Satellite(row)
            self.constellation[sat.get_uuid] = sat

        
   # Creates a memmap to read parsed satellite data. Note that this function assumes the parser
   # has already been ran by itself.
    def load_data_matrix(self, satName):
        filename = 'parsed_' + satName + '.npy'
        if filename in os.listdir(self.parsed_datapath):
            matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
        else:
            matrix = 0
            
        return np.array(matrix)

        
    def load_sensor_model(self, satName):
        filename = 'sensor_parsed_' + satName + '.json'
        with open(parsed_datapath + filename, "w") as sensor_json:
                sensorModel = json.load(sensor_json)       
        return sensorModel

        
    def load_extra_info(self):
        ei_Matrices = []

        for filename in os.listdir(self.parsed_datapath):
            if(filename[0] == 'e'):
                matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
                ei_Matrices.append(matrix)
            
        return np.array(ei_Matrices)    