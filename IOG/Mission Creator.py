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
import numpy as np

import time

class MissionCreator:


    #initialization function for the mission creator class    
    def __init__(self, datapath, missionpath, parsed_datapath):
        self.test = False

        t0 = time.time()
        self.datapath = datapath
        self.missionpath = missionpath
        self.parsed_datapath = parsed_datapath
        
        self.system_setup(datapath, missionpath) 
        
        missions = self.parser.create_missions()
        dataMatrices = self.load_data_matrices()
        extraInfoMatrix = self.load_extra_info()
        
        #Data Loading Test Code
        if(self.test):
            print(len(dataMatrices))
            print(len(dataMatrices[0]))
            print(len(dataMatrices[0][0]))
        
        for mission in missions:        
           self.generate_imaging_opportunities(mission, dataMatrices, extraInfoMatrix) #main functions
        t1 = time.time()
        deltaT = t1-t0
        print('Total Time: ' + str(deltaT))
        
        #Mission test code
        if(self.test):        
            print(str(len(missions)))
            print(missions[0].get_name())
            print(missions[0].get_coordinates())
            print(missions[0].get_interval_start_time())

    
    # Initializes blocks of one hierarchy level lower that this one.
    # At the moment, the Constellation is unused.
    def system_setup(self, datapath, misionpath):
        self.parser = Parser(self.datapath, self.missionpath)
        self.constellation = Constellation()
        self.calculator = Calculator()
    

    # Starts the calculator's processes to generate results    
    def generate_imaging_opportunities(self, mission, dataMatrices, extraInfoMatrix):
        self.calculator.generate_imaging_opportunities(mission, dataMatrices, extraInfoMatrix)


    # Creates a memmap to read parsed satellite data. Note that this function assumes the parser
    # has already been ran by itself.
    # This process will eventually be moved to the Satellite/Constellation classes 
    def load_data_matrices(self):
        dataMatrices = []

        for filename in os.listdir(self.parsed_datapath):

            if(filename[0] == 'p'):
                matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
                dataMatrices.append(matrix)
            
        return np.array(dataMatrices)
        
        
    def load_extra_info(self):
        ei_Matrices = []

        for filename in os.listdir(self.parsed_datapath):

            if(filename[0] == 'e'):
                matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')
                ei_Matrices.append(matrix)
            
        return np.array(ei_Matrices)    
    
#main code goes here
#we will eventually define a main function so that the program can run as an executable

#The two parameters are the filepaths for the satellite data and the mission info data, respectively
if __name__ == "__main__":
    missionCreator = MissionCreator("../../Data/", "../../Missions/Alpha Test/", "../../Parsed Data/")

