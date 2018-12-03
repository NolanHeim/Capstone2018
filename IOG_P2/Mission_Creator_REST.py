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
from Constellation import *
from Mission import *
import numpy as np
import json
import uuid

import time

class Mission_Creator_REST:


    #initialization function for the mission creator class    
    def __init__(self, parsed_datapath, results_path):
        self.test = False
        self.results_path = results_path
        self.parsed_datapath = parsed_datapath
        self.calculator = Calculator()
        self.constellation = Constellation(self.parsed_datapath)
    

    # Starts the calculator's processes to generate results    
    def generate_imaging_opportunities(self, input_json, uuid):

        t0 = time.time()
        print("in MC")

        mission = self.create_mission_from_json(input_json, uuid)
        print("mission made")

        ## Constellation add all mission satellite UUID's not already in.
        if (mission.check_params() == False):
            return "ERR"

        print("got passed check params")
        uuid_list = mission.get_ids_to_consider()
        self.constellation.add_satellite_data(uuid_list)
        
        print("satellites added")
        #this should be a list of a list of windows: where the first list corresponds to satellites,
        #and the list within it to the windows
        [windows_per_sat, sats] = self.calculator.generate_imaging_opportunities(mission, self.constellation)

        print("out of calc")

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
        input_dict = input_json
        targetCoordinates = input_dict.get("TargetArea", "")
        print("got area")
        name = str(uuid)
        startTime = input_dict.get("POI", "").get("startTime", "")
        endTime = input_dict.get("POI", "").get("endTime", "")
        
        print("got POI")
        if("SensorID" in input_dict):
            idsToConsider = input_dict.get("PlatformID")
        else:
            idsToConsider = []
        
        sensorType = ""
        minSolarAngle = 0.0
        maxSolarAngle = 180.0
        
        if("Filter" in input_dict):
            if("SolarAngles" in input_dict.get("Filter")):
                minSolarAngle = input_dict.get("Filter").get("SolarAngles").get("MinimumIncidenceAngle")
                maxSolarAngle = input_dict.get("Filter").get("SolarAngles").get("MaximumIncidenceAngle")
            if("SensorType" in input_dict.get("Filter")):
                sensorType = input_dict.get("Filter").get("SensorType")
            
                    
        print("got filter")
        mission = Mission(targetCoordinates, name, sensorType, startTime, endTime, 
                          idsToConsider, minSolarAngle, maxSolarAngle)
        return mission        


if __name__ == '__main__':
    with open("test_mission.json", "r") as missionjson:
        missiondict = json.load(missionjson)
    
    #missioninput = json.dumps(missiondict)
    newid = uuid.uuid4()
    
    MC = Mission_Creator_REST("../../Parsed Data/", "../../Saved Results")
    MC.generate_imaging_opportunities(missiondict, newid)
