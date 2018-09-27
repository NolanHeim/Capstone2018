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

###
# Example Call:
#X = Parser('../../Data')
#X.parseData('RCM1_Reference_Orbit.e')
#Y = X.getDataMatrix()
###

class Parser:
    
    def __init__(self, path):
        self.path = path;
    
    def parseData(self, file):
        dataFile = os.path.join(self.path, file)
        fileObj = open(dataFile, 'r')
        fileLines = [line.rstrip('\n') for line in fileObj]
        self.lines = fileLines
        fileObj.close()
        
        dataStartIndex = fileLines.index('EphemerisTimePosVel') + 1
        dataStopIndex = fileLines.index('END Ephemeris') - 1
        self.dataMatrix = [[float(i) for i in fileLines[j].split()] for j in range(dataStartIndex,dataStopIndex)]
        
    def getDataMatrix(self):
        return self.dataMatrix



