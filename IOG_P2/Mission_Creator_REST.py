# -*- coding: utf-8 -*-
#
# Mission Creator.py
#
# This is the highest level class in terms of project hierarchy.
# This class initializes other classes as well as reads mission parameters.
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#

from Parser import *
from Calculator import *
import numpy as np
import json

import time

class Mission_Creator_REST:


    #initialization function for the mission creator class    
    def __init__(self, parsed_datapath, results_path):
        self.test = False
        self.results_path = results_path
        self.parsed_datapath = parsed_datapath
        self.calculator = Calculator()
    

    # Starts the calculator's processes to generate results    
    def generate_imaging_opportunities(self, input_json, uuid):
        dataMatrices = self.load_data_matrices()
        extraInfoMatrix = self.load_extra_info()

        t0 = time.time()

        mission = self.create_mission_from_json(input_json, uuid)
        
        if (mission.check_params() == False):
            return "ERR"
        
        #mission.display_parameters()
        
        #this should be a list of a list of windows: where the first list corresponds to satellites,
        #and the list within it to the windows
        [windows_per_sat, sats] = self.calculator.generate_imaging_opportunities(mission, dataMatrices, extraInfoMatrix)

        #print(len(windows_per_sat))
        #print(sats)
        #print(len(windows_per_sat[0]))
        #print(windows_per_sat[0])

        opportunity_jsons = []
        for i in range(0, len(sats)):
            #create an opportunity json for this satellite
            opportunities = {
                "Platform ID": str(sats[i]),
                "Oppotunity Window": str(windows_per_sat[i]),
            }

            
            opp_json = json.dumps(opportunities)
            opportunity_jsons.append(opp_json)

        t1 = time.time()
        deltaT = t1-t0
        print('Total Time: ' + str(deltaT))
        
        #print(opportunity_jsons)        
        
        return opportunity_jsons


    #makes a mission object based on the queried visibility search input json
    def create_mission_from_json(self, input_json, uuid):
        #print(json.dumps(input_json))
        #print("making dictionary")
        input_dict = input_json
        print(input_dict)
        targetCoordinates = input_dict.get("Target", "")
        print(targetCoordinates)
        name = str(uuid)
        startTime = input_dict.get("POI", "").get("startTime", "")
        endTime = input_dict.get("POI", "").get("endTime", "")
        
        if("PlatformID" in input_dict):
            idsToConsider = input_dict.get("PlatformID")
        else:
            idsToConsider = []
            
        mission = Mission(targetCoordinates, name, "", "", 0, startTime, endTime, idsToConsider)
        return mission        


    # Creates a memmap to read parsed satellite data. Note that this function assumes the parser
    # has already been ran by itself.
    # This process will eventually be moved to the Satellite/Constellation classes 
    def load_data_matrices(self):
        dataMatrices = []

        for filename in os.listdir(self.parsed_datapath):
            if(filename[0] == 'p'):
                matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
                dataMatrices.append(matrix)
            
        return np.array(dataMatrices)
        
        
    def load_extra_info(self):
        ei_Matrices = []

        for filename in os.listdir(self.parsed_datapath):

            if(filename[0] == 'e'):
                matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
                ei_Matrices.append(matrix)
            
        return np.array(ei_Matrices)    
    