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

class MissionCreator:


    #initialization function for the mission creator class    
    def __init__(self, datapath, missionpath, parsed_datapath):

        self.datapath = datapath
        self.missionpath = missionpath
        self.parsed_datapath = parsed_datapath
        #initialize dependent objects
        self.system_setup(datapath, missionpath) 
        
        #parser already contains missionpath
        missions = self.parser.create_missions()
        #dataMatricesWithExtraInfo = np.array(self.load_data_matrices())
        
        #dataMatrix = []
        #Extract the top vector of each data matrix.
        #for dataline in dataMatricesWithExtraInfo[:,0]:
        #    dataMatrix.append(dataline)
        
        dataMatrices = self.load_data_matrices()

        #Delete the extra top vector from the data matrices
        # TODO Might need to be modified for multiple files.
        #dataMatrices = np.array(dataMatricesWithExtraInfo[0,1])
        
        #Data Loading Test Code
        #print(len(dataMatrices))
        #print(len(dataMatrices[0]))
        #print(len(dataMatrices[0][0]))
        
        #print("STOP")        
        
        for mission in missions:        
           self.generate_imaging_opportunities(mission, dataMatrices) #main functions
        
        #Mission test code        
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

    
    def load_data_matrices(self):
        dataMatrices = []
        #self.parser.parse_data(self.parsed_datapath)
        #print("MC PDP is "+self.parsed_datapath)
        for filename in os.listdir(self.parsed_datapath):
            #dataMatrices.append(self.load(self.parsed_datapath+filename))
            
            #read_map = np.memmap(filename, dtype='float32', mode='w+', shape=)
            matrix = np.load(self.parsed_datapath+filename, mmap_mode='r')

            dataMatrices.append(matrix)

#            dataMatrices.append(np.load(self.parsed_datapath+filename))
            #dataMatrices = (np.load(self.parsed_datapath+filename))
            
        return np.array(dataMatrices)
    
    
#main code goes here
#we will eventually define a main function so that the program can run as an executable

#The two parameters are the filepaths for the satellite data and the mission info data, respectively
missionCreator = MissionCreator("../../Data/", "../../Missions/", "../../Parsed Data/")

