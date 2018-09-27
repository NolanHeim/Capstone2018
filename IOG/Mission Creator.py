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

import numpy
import matplotlib
import requests
import flask

class MissionCreator:
    calculator;
    constellation;
    #TODO insert class functions and parameters

    #initialization function for the mission creator class    
    def __init__(self, datapath, missionpath):
        print(datapath)
        print(missionpath)
        
        systemSetup(datapath)        
        
        while(missionsExist(missionpath)):
            currentMission = loadCurrentMission()
            generateImagingOpportunities(currentMission)
        #initialize variables
        
    def systemSetup(datapath):
        self.calculator = Calculator(datapath);
        self.constellation = Constellation(datapath);
        
    def loadCurrentMission():
        
        
    def missionsExist(datapath):
        
    
    def generateImagingOpportunities(mission):
        calculator.generateImagineOpportunities(mission, constellation)
        

#main code goes here
#we will eventually define a main function so that the program can run as an executable

#The two parameters are the filepaths for the satellite data and the mission info data, respectively
missionCreator = MissionCreator("this very folder", "this folder's parent")

