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
# 2. Parses sensor data files (.txt) with sensor data such as positional offset
#    and angular freedom in each direction, as well as sensor type and whether
#    or not sensor angles are dependent.
# 
# Inputs:
#   orbitdatapath: The relative location of the satellite (orbit) data. (ex: ../../Data/Orbit/)
#   sensordatapath: Relative location of the sensor data. (ex: ../../Data/Sensor/)
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


class Parser:
    
    def __init__(self, orbitdatapath, sensordatapath):
        self.orbitdatapath = orbitdatapath        
        self.sensordatapath = sensordatapath
        self.jDate_offset = 3
        self.SENSOR_INDEX_INCREMENT = 6
    
    
   # This function parses all required data from satellite .e files.
   # It outputs the orbital data (times, positions, velocities) and extra info
   # (JDate, UUID) to the disk, in the folder specified by parsed_datapath
    def parse_orbit_data(self, parsed_datapath):
        
        for filename in os.listdir(self.orbitdatapath): 
            
            # Open the file, parse the contents
            dataFile = os.path.join(self.orbitdatapath, filename)
            fileObj = open(dataFile, 'r')
            fileLines = [line.rstrip('\n').strip(' ').strip('\t') for line in fileObj]
            fileLines[:] = [x for x in fileLines if x != '']
            fileObj.close()

            # Find the JDate in the file        
            JDateIndex = fileLines.index('BEGIN Ephemeris') + self.jDate_offset
            JDLine = fileLines[JDateIndex].split()
            jDate = float(JDLine[-1])
            
            # Generate a UUID for this satellite
            satelliteID = str(uuid.uuid4())        
            satelliteName = os.path.splitext(filename)[0]
            
            # Find the start of the orbital data (times, positions, velocities)
            # and generate the corresponding matrix
            dataStartIndex = fileLines.index('EphemerisTimePosVel') + 1
            dataStopIndex = fileLines.index('END Ephemeris') - 1
            dataMatrix = [[float(i) for i in fileLines[j].split()] for j in range(dataStartIndex,dataStopIndex)]
            
            # Paths and filenames for the orbital and extra info files that will be created.
            parsed_dataFile = os.path.join(parsed_datapath, 'parsed_'+satelliteName)
            extraInfoFileName = os.path.join(parsed_datapath, 'ex_parsed_'+satelliteName)
            
            # Implement the extra info as a dictionary, then save as json for easy
            # extraction later. Saving as a numpy array like in Phase 1 did not work
            # due to conflicting data types.
            extraInfo = {}
            extraInfo["jDate"] = jDate
            extraInfo["satelliteID"] = satelliteID
            extraInfo["satelliteName"] = satelliteName
            extraInfo["parsed_dataFile"] = 'parsed_'+satelliteName+'.npy'

            # save the extra info
            with open(extraInfoFileName+".json", "w") as extraInfoFile:
                extraInfoFile.write(json.dumps(extraInfo))
        
            # Save the data matrix
            data = np.array(dataMatrix)
            np.save(parsed_dataFile, data)
            

    # parses the sensor data for each satellite
    # 
    # Inputs:
    #   parsed_datapath - path for satellite data to be saved at.
    def parse_sensor_data(self, parsed_datapath):
        
        # There should exist no more than 1 file per satellite, with every sensor
        # for a given satellite contained in the file. Sensor file structure should
        # also follow the specification given in the final report.
        for filename in os.listdir(self.sensordatapath):
            sensorDict = {}
            sensorDataFile = os.path.join(self.sensordatapath, filename)
            fileObj = open(sensorDataFile, 'r')
            fileLines = [line.rstrip('\n') for line in fileObj]
            fileLines[:] = [x for x in fileLines if x != '']
            fileObj.close()
            satelliteName = os.path.splitext(filename)[0]

            #Figures out how many sensors there are on this satellite
            numSensors = int(fileLines[0].split(" ")[3])
            baseI = 1

            # Repeats for each sensor in the file for this satellite
            for sensorIndex in range(1,numSensors+1):
                sensorName = ""
                nameLine = fileLines[baseI].split(" ")
                for index in range(1, len(nameLine)):
                    sensorName = sensorName + nameLine[index] + " "
                sensorName = sensorName.strip(" ")
                
                sensorDict[sensorName] = {}
                sensorDict[sensorName]['Type'] = fileLines[baseI+1].split(' ')[1]
                sensorDict[sensorName]['Imaging Type'] = fileLines[baseI+2].split(' ')[2]

                # Find and store orientation angles for this sensor                
                orientation = [float(i) for i in fileLines[baseI+3].split(' ')[1].split(',')]
                
                sensorDict[sensorName]['Orientation'] = orientation
                sensorDict[sensorName]['Angle Dependent'] = fileLines[baseI+4].split(' ')[2]

                if(fileLines[baseI+4].split(' ')[2] == 'True'):
                    #Expect 1 Angles
                    sensorDict[sensorName]['Rotation'] = float(fileLines[baseI+5].split(' ')[1])                    
                else:
                    #Expect 2 Angles
                    rotation = [float(i) for i in fileLines[baseI+5].split(' ')[1].split(',')]
                    sensorDict[sensorName]['Rotation'] = rotation
                
                #Increment baseI for the next sensor
                baseI = baseI + self.SENSOR_INDEX_INCREMENT
                
                #NOTE: Could add a check here to see if they are trying a dependent push broom (can't happen)
                # AND/OR: Could check to see if the expected 1 or 2 angles is incorrect. 
            with open(parsed_datapath + "sensor_parsed_"+satelliteName+".json", "w") as sensor_output:
                sensor_output.write(json.dumps(sensorDict))
        
        
    # Top level function, controls the parsing of Satellite and Sensor data
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
            

# Code that runs the parser as standalone. This just parses and stores satellite and sensor data. 
if __name__ == "__main__":
    parser = Parser('../../Data/Orbit/', '../../Data/Sensor/')
    parser.parse_data('../../Parsed Data/')

