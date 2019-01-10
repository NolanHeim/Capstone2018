# -*- coding: utf-8 -*-
#
# Mission_Creator_REST.py
#
# (PHASE 2)
#
# This is the highest level class in terms of project hierarchy.
# This class initializes other classes as well as reads mission parameters.
#
# Note that it is called Mission_Creator_REST to differentiate from the legacy
# file Mission Creator.py, which was written before the REST API was implemented
# and is no longer in use.
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#

from Parser import *
from Calculator import *
from Constellation import *
from Mission import *
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

        mission = self.create_mission_from_json(input_json, uuid)

        if (mission.check_params() == False):
            return "ERR"

        uuid_list = mission.get_ids_to_consider()

        ## Constellation add all mission satellite UUID's not already in.
        self.constellation.add_satellite_data(uuid_list)
        
        # windows_per_sat is a list of a list of windows: where the first list corresponds to 
        # the windows for each satellite, and the list within it to the actual windows.
        # sats is the corresponding list of satellites.
        [windows_per_sat, sats] = self.calculator.generate_imaging_opportunities(mission, self.constellation)


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
        
        
        return opportunity_jsons


    #makes a mission object based on the queried visibility search input json
    def create_mission_from_json(self, input_json, uuid):
        input_dict = input_json
        targetCoordinates = input_dict.get("TargetArea", "")
        name = str(uuid)
        startTime = input_dict.get("POI", "").get("startTime", "")
        endTime = input_dict.get("POI", "").get("endTime", "")
        
        #print("got POI")
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
            
                    
        mission = Mission(targetCoordinates, name, sensorType, startTime, endTime, 
                          idsToConsider, minSolarAngle, maxSolarAngle)
        return mission        


# For debugging in-IDE only. This file is usually only called through API.py.
if __name__ == '__main__':
    with open("test_mission.json", "r") as missionjson:
        missiondict = json.load(missionjson)
    
    newid = uuid.uuid4()
    
    MC = Mission_Creator_REST("../../Parsed Data/", "../../Saved Results")
    MC.generate_imaging_opportunities(missiondict, newid)
