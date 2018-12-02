# -*- coding: utf-8 -*-
#
# Parser.py
#
# Parser has two main functions:
#
# 1. Parses "ephemeris" (.e) files with satellite data. This occurs when this file is run
#  standalone. Saves a numpy array in the 
#  Assumes the files have "EphemerisTimePosVel" before the satellite orbital data
#   and "END Ephemeris" at the end of the satellite data.
#  Assumes data is in the format:
#   Time (s), X (m), Y(m), Z(m), Vx(m), Vy(m), Vz(m)
#   Where X, Y, and Z are the positions and Vx, Vy and Vz are velocities.
#
# 2. Parses mission files (.txt) stored in the directory specified by missionpath. 
#  Creates corresponding mission objects that are returned as a list. This functionality
#  is not ran when the code is run standalone. The function to do this gets called within
#  the Mission Creator class.
# 
# Inputs:
#   orbitdatapath: The relative location of the satellite (orbit) data. (ex: ../../Data/Orbit/)
#   sensordatapath: Relative location of the sensor data. (ex: ../../Data/Sensor/)
#   missionpath: The relative location of all satellite mission files (ex: ../../Missions/)
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
import os
import numpy as np
from Mission import *
import time
import uuid
import json

###IMPORTANT - CHANGE THIS WHEN WE MOVE ON FROM PHASE 1###
phase1 = True
###


class Parser:
    
    def __init__(self, orbitdatapath, sensordatapath):
        self.orbitdatapath = orbitdatapath        
        self.sensordatapath = sensordatapath
        
        #self.lineLeaders = ['Name:', 'Sensor_Type:', 'Illumination_Direction:',
        #                    'Illumination_Threshold:', 'Interval_Start_Time:',
        #                    'Interval_End_Time:']
        #self.p_name = 0
        #self.p_sensorType = 1
        #self.p_illumDir = 2
        #self.p_illumThresh = 3
        #self.p_intervalStart = 4
        #self.p_intervalEnd = 5
        
        #self.jDate_offset = 3
        
        #self.INDEX_OF_COORDINATES = 3
        #self.SENSOR_INDEX_INCREMENT = 5
        #self.databaseFilename = "Satellite_Database.json"
    
    
    # Parses the orbital data for each satellite
    #
    # Inputs:
    #   parsed_datapath - path for satellite data to be saved at.
    #
    def parse_orbit_data(self, parsed_datapath):
        for filename in os.listdir(self.orbitdatapath): 
            dataFile = os.path.join(self.orbitdatapath, filename)
            fileObj = open(dataFile, 'r')
            fileLines = [line.rstrip('\n') for line in fileObj]
            fileLines[:] = [x for x in fileLines if x != '']
            fileObj.close()
        
            JDateIndex = fileLines.index('BEGIN Ephemeris') + self.jDate_offset
            JDLine = fileLines[JDateIndex].split()
            jDate = float(JDLine[-1])
            
            satelliteID = uuid.uuid4()        
            satelliteName = os.path.splitext(filename)[0]
            self.updateDatabase(satelliteName, satelliteID, parsed_datapath)
        
            dataStartIndex = fileLines.index('EphemerisTimePosVel') + 1
            dataStopIndex = fileLines.index('END Ephemeris') - 1
            dataMatrix = [[float(i) for i in fileLines[j].split()] for j in range(dataStartIndex,dataStopIndex)]
            
            newFile = os.path.join(parsed_datapath, 'parsed_'+satelliteName)
            extraInfoFile = os.path.join(parsed_datapath, 'ex_parsed_'+satelliteName)
            
            extraInfo = [jDate, satelliteID, satelliteName, newFile, sensorFile]            
            
            data = np.array(dataMatrix)
            np.save(newFile, data)
            
            extra = np.array(extraInfo)
            np.save(extraInfoFile, extra)


    # parses the sensor data for each satellite
    # 
    # Inputs:
    #   parsed_datapath - path for satellite data to be saved at.
    #
    def parse_sensor_data(self, parsed_datapath):
        
        for filename in os.listdir(self.sensordatapath):
            sensorDict = {}
            sensorDataFile = os.path.join(self.sensordatapath, filename)
            fileObj = open(sensorDataFile, 'r')
            fileLines = [line.rstrip('\n') for line in fileObj]
            fileLines[:] = [x for x in fileLines if x != '']
            fileObj.close()
            satelliteName = os.path.splitext(filename)[0]
            numSensors = fileLines[0][3] 
            baseI = 1
            for sensorIndex in range(1,numSensors+1):
                sensorName = fileLines[baseI][1:]
                sensorDict[sensorName] = {}
                sensorDict[sensorName]['Type'] = fileLines[baseI+1][1:]
                orientation = [float(i) for i in fileLines[baseI+2][1:].split(',')]
                sensorDict[sensorName]['Orientation'] = orientation
                sensorDict[sensorName]['Angle Dependent'] = fileLines[baseI+3][2]
                if(fileLines[baseI+3][2] == 'True'):
                    #Expect 1 Angles
                    sensorDict[sensorName]['Rotation'] = float(fileLines[baseI+4][1])                    
                else:
                    #Expect 2 Angles
                    rotation = [float(i) for i in fileLines[baseI+4][1:].split(',')]
                    sensorDict[sensorName]['Rotation'] = rotation
                #Increment baseI for the next sensor
                baseI = baseI + self.SENSOR_INDEX_INCREMENT
                #NOTE: Could add a check here to see if they are trying a dependent push broom (can't happen)
                # AND/OR: Could check to see if the expected 1 or 2 angles is incorrect. 
            with open(parsed_datapath + "sensor_parsed_"+satelliteName+".json", "w") as sensor_output:
                sensor_output.write(json.dumps(sensorDict))
        
    # Controls the Parsing of Satellite data
    #
    # Inputs:
    #   parsed_datapath - path for satellite data to be saved at.
    #
    def parse_data(self, parsed_datapath):
        t0 = time.time()
        #Parse orbital data & save
        self.parse_orbit_data(parsed_datapath)
        #Parse sensor data & save2
        self.parse_sensor_data(parsed_datapath)
        print("Total time = "+str(time.time() - t0))
            

    # Parses mission files stored as .txt files in folder specified by missionpath.
    # Note that to work correctly, mission files must be compliant with the guidelines outlined
    # in the Mission Specifications. 
    #def create_missions(self):
    #    missions = []        
        
    #    for filename in os.listdir(self.missionpath):
    #        fileObj = open(self.missionpath + filename, 'r')
    #        fileLines = [line.rstrip('\n') for line in fileObj]
    #        fileObj.close()
            
    #        targetCoordinates = []
    #        for i in range(1,self.coordinateIndexCap+1):
    #            targetCoordinates.append(self.get_mission_coordinates(fileLines,i))
    #    
    #        name = self.get_simple_param(fileLines, self.p_name)
    #        sensorType = self.get_simple_param(fileLines, self.p_sensorType).lower()
    #        illumDir = self.get_simple_param(fileLines, self.p_illumDir).lower()
    #        illumThresh = float(self.get_simple_param(fileLines, self.p_illumThresh))
    #        intervalStart = self.get_simple_param(fileLines, self.p_intervalStart)
    #        intervalEnd = self.get_simple_param(fileLines, self.p_intervalEnd)            
    #        idsToConsider = []
    #        
    #        newMission = Mission(targetCoordinates, name, sensorType, illumDir, illumThresh, intervalStart, intervalEnd, idsToConsider)
    #        
    #        missions.append(newMission)
    #        
    #    return missions
        
        
    # Returns the coordinates of the mission as an ordered pair representing [NS, EW] coordinates. For 
    # phase 1, only one coordinate pair is returned (the target is a point). For phase 2, this will return
    # four coordinates which will form a target window.
    #def get_mission_coordinates(self, fileLines, index):
    #    #first, get the coordinates in string format
    #    
    #    for i in range(0, len(fileLines)):
    #        if(len(fileLines[i]) != 0 and fileLines[i].split()[0] == 'Coordinate' and fileLines[i].split()[2] == str(index)+':'):
    #            coordinateStringLine = fileLines[i].split()
        
    #    coordinateString = coordinateStringLine[self.INDEX_OF_COORDINATES]
        
    #    #next, convert the coordinate into actual numbers
    #    coordinateStringCommaIndex = coordinateString.index(',')        
    #    
    #    coordinateStringNS = coordinateString[2:coordinateStringCommaIndex] 
    #    coordinateStringEW = coordinateString[coordinateStringCommaIndex+2:-1]
    #    
    #    coordinateNS = float(coordinateStringNS)
    #    coordinateEW = float(coordinateStringEW)
    #            
    #    if coordinateString[1] == '-':
    #        coordinateNS = -coordinateNS
    #        
    #    if coordinateString[coordinateStringCommaIndex + 1] == '-':
    #        coordinateEW = -coordinateEW
    #    
    #    coordinates = []
    #    coordinates.append(coordinateNS)
    #    coordinates.append(coordinateEW)
    #    
    #    return coordinates


    # Gets any param other than coordinates, eg name, date, etc. 
    #def get_simple_param(self, fileLines, param):
    #    line = ''

    #    for i in range(0, len(fileLines)):
    #        if(len(fileLines[i]) != 0 and fileLines[i].split()[0] == self.lineLeaders[param]):
    #            line = fileLines[i].split()        
    #    
    #    if(len(line) != 0):
    #        if(param == self.p_intervalStart or param == self.p_intervalEnd):
    #            return line[1]+' '+line[2]
    #        else:
    #            return line[1]
    #    else:
    #        return 'ERROR'


    #keeps track of which satellites we have
    #
    # Don't think this is currently working...
    #def updateDatabase(self, satelliteName, satelliteID, parsed_datapath):
    #    if os.path.isfile(parsed_datapath + self.databaseFilename): #then our database already exists
    #        with open(parsed_datapath + self.databaseFilename, "w") as database_json:
    #            database_dict = json.load(database_json)
    #            if(satelliteID not in database_dict):
    #                database_dict.update({satelliteID : satelliteName})
    #            json.dump(database_dict, database_json)
    #    else:
    #        with open(parsed_datapath + self.databaseFilename, "w") as database_json:
    #            new_dict = {satelliteID : satelliteName}
    #            json.dump(database_json, new_dict)
                

# Code that runs the parser as standalone. This just parses and stores satellite data. 
if __name__ == "__main__":
    parser = Parser('../../Data/', '../../Missions/Alpha Test/')
    parser.parse_data('../../Parsed Data/')

