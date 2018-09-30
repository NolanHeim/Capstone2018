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
from Constellation import *
from Calculator import *


class MissionCreator:


    #initialization function for the mission creator class    
    def __init__(self, datapath, missionpath):

        self.datapath = datapath
        self.missionpath = missionpath
        #initialize dependent objects
        self.system_setup(datapath, missionpath) 
        
        #parser already contains missionpath
        missions = self.parser.create_missions()
        dataMatrices = self.parser.make_data_matrices()       
        
        for mission in missions:        
            self.generate_imaging_opportunities(mission, dataMatrices) #main functions
        print(str(len(missions)))
        print(missions[0].get_name())
        print(missions[0].get_coordinates())
        print(missions[0].get_interval_start_time())
    
    #initializes blocks of one hierarchy level lower that this one   
    def system_setup(self, datapath, misionpath):
        self.parser = Parser(self.datapath, self.missionpath)
        self.constellation = Constellation()
        self.calculator = Calculator()
    
    
    def generate_imaging_opportunities(self, mission, dataMatrices):
        self.calculator.generate_imaging_opportunities(mission, dataMatrices)
        print('hi')

#main code goes here
#we will eventually define a main function so that the program can run as an executable

#The two parameters are the filepaths for the satellite data and the mission info data, respectively
missionCreator = MissionCreator("../../Data", "../../Missions/")

