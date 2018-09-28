# -*- coding: utf-8 -*-
#
# Parser.py
#
# Parses "ephemeris" (.e) files with satellite data.
# Assumes the files have "EphemerisTimePosVel" before the satellite orbital data
#   and "END Ephemeris" at the end of the satellite data.
# Assumes data is in the format:
#   Time (s), X (m), Y(m), Z(m), Vx(m), Vy(m), Vz(m)
#   Where X, Y, and Z are the positions and Vx, Vy and Vz are velocities.
#
# Inputs:
#   path: The relative location of the satellite data. (ex: ../../Data)
#   file: The file name. (ex: 'RCM1_Reference_Orbit.e')
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
import os

###IMPORTANT - CHANGE THIS WHEN WE MOVE ON FROM PHASE 1###
phase1 = True
###

###
# Example Call:
#X = Parser('../../Data')
#X.parse_data('RCM1_Reference_Orbit.e')
###

class Parser:
    
    def __init__(self, datapath, missionpath):
        self.datapath = datapath
        self.missionpath = missionpath
    
        #how many coordinates there are
        if(phase1):
            self.coordinateIndexCap = 1
        else:
            self.coordinateIndexCap = 4        
            
        self.INDEX_OF_COORDINATES = 4
        self.INDEX_OF_NAME = 2
        self.INDEX_OF_SENSOR_TYPE = 3
        
        self.lineLeaders = ['Name:', 'Sensor_Type:', 'Illumination_Direction:',
                            'Illumination_Threshold:', 'Interval_Start_Time:',
                            'Interval_End_Time:']
        self.p_name = 0
        self.p_sensorType = 1
        self.p_illumDir = 2
        self.p_illumThresh = 3
        self.p_intervalStart = 4
        self.p_intervalEnd = 5
   
    def parse_data(self, file):
        dataFile = os.path.join(self.path, file)
        fileObj = open(dataFile, 'r')
        fileLines = [line.rstrip('\n') for line in fileObj]
        self.lines = fileLines
        fileObj.close()
        
        dataStartIndex = fileLines.index('EphemerisTimePosVel') + 1
        dataStopIndex = fileLines.index('END Ephemeris') - 1
        dataMatrix = [[float(i) for i in fileLines[j].split()] for j in range(dataStartIndex,dataStopIndex)]
        
        return dataMatrix


    def create_missions(self):
        missions = []        
        
        for filename in os.listdir(self.missionpath):
            fileObj = open(filename, 'r')
            fileLines = [line.rstrip('\n') for line in fileObj]
            fileObj.close()
            
            targetCoordinates = []
            for i in range(1,self.coordinateIndexCap):
                targetCoordinates.append(self.get_mission_coordinates(fileLines,i))
        
            name = self.get_simple_param(fileLines, self.p_name)
            sensorType = get_simple_param(fileLines, self.p_sensorType).lower()
            illumDir = self.get_simple_param(fileLines, self.p_illumDir).lower()
            illumThresh = float(self.get_simple_param(fileLines, self.p_illumThresh))
            intervalStart = self.get_simple_param(fileLines, self.p_intervalStart)
            intervalEnd = self.get_simple_param(fileLines, self.p_intervalEnd)            
            
            newMission = Mission(targetCoordinates, name, sensorType, illumDir, illumThresh, intervalStart, intervalEnd)
            
            missions.append(newMission)
            
        return missions
        
        
    def get_mission_coordinates(self, fileLines, index):
        #first, get the coordinates in string format
        coordinateStringIndex = fileLines.index('Coordinate Corner '+str(index)+':') #gets the line that we care about
        coordinateStringLine = fileLines[coordinateStringIndex].split() #splits  based on whitespaces
        coordinateString = coordinateStringLine[self.INDEX_OF_COORDINATES]
        
        #next, convert the coordinate into actual numbers
        coordinateStringCommaIndex = coordinateString.index(',')        
        
        coordinateStringNS = coordinateString[2:coordinateStringCommaIndex - 1] 
        coordinateStringEW = coordinateString[coordinateStringCommaIndex+1:-1]
        
        coordinateNS = float(coordinateStringNS)
        coordinateEW = float(coordinateStringEW)
                
        if coordinateString[1] == '-':
            coordinateNS = -coordinateNS
            
        if coordinateString[coordinateStringCommaIndex + 1] == '-':
            coordinateEW = -coordinateEW
        
        coordinates = []
        coordinates.append(coordinateNS)
        coordinates.append(coordinateEW)
        
        return coordinates


    def get_simple_param(self, fileLines, param):
        lineIndex = fileLines.index(self.lineLeaders(param))
        line = fileLines[lineIndex].split()
        
        return line[1]
