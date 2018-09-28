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


    #initialization function for the mission creator class    
    def __init__(self, datapath, missionpath):
        print(datapath) #for testing, remove after
        print(missionpath) #for testing, remove after
        
        #initialize dependent objects
        system_setup(datapath, missionpath) 
        
        #parser already contains missionpath
        missions = parser.make_missions()
        for mission in missions:        
            generate_imaging_opportunities(mission) #main functions
        
    
    #initializes blocks of one hierarchy level lower that this one   
    def system_setup(datapath, misionpath):
        self.parser = Parser(datapath, missionpath);
        self.constellation = Constellation(self.parser);
        self.calculator = Calculator(self.constellation);
        
    #passes params into calculator for one mission
    def load_current_mission():
        
        
    #sees if there are still missions for which outputs have not been calculated
    def missions_exist(datapath):
        return 1
        
    
    def generate_imaging_opportunities(mission):
        calculator.generate_imaging_opportunities(mission)
        

#main code goes here
#we will eventually define a main function so that the program can run as an executable

#The two parameters are the filepaths for the satellite data and the mission info data, respectively
missionCreator = MissionCreator("this very folder", "this folder's parent")

