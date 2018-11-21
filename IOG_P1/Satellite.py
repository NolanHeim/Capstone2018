# -*- coding: utf-8 -*-
#
# Satellite.py
#
# Represents the data for a single satellite
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
from Parser import *
from bisect import bisect_left
import math

class Satellite:
    #TODO insert class functions and parameters    
    def __init__(self, name, path, dataFile):
        self.dataFile = dataFile
        self.name = name
        dataParser = Parser(path)
        dataParser.parseData(dataFile)
        self.dataMatrix = dataParser.getDataMatrix()
        
    def binaryListSearch(self,dataList, target):
        index = bisect_left(dataList, target)
        if index == 0:
            return 0
        if index == len(dataList):
            return len(dataList)
        beforeIndex = dataList[index - 1]
        afterIndex = dataList[index]
        if afterIndex - target < target - beforeIndex:
            return index-1
        else:
            return index        
        
        
    #Assumes time is in seconds
    def getRadius(self,time):
        timeVector = [row[0] for row in self.dataMatrix]
        rowIndex = self.binaryListSearch(timeVector, time)
        radius = math.sqrt(self.dataMatrix[rowIndex][1]**2 + self.dataMatrix[rowIndex][2]**2 + self.dataMatrix[rowIndex][3]**2)
        return radius
        
    def getSatelliteData(self):
        return self.dataMatrix