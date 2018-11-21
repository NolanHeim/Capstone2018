# -*- coding: utf-8 -*-
#
# Constellation.py
#
# Represents the group of all satellites for which data is provided.
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
from Satellite import *

class Constellation:    
    
    #TODO insert class functions and parameters
    def __init__(self):
        self.constellation = [];
        self.satellites = [];
        #initialize constellations
        
        
    #gets the info for one satellite
    def add_to_constellation(self, satellite):
        self.constellation.append(satellite)

    def new_satellite(self, name, dataFile):
        self.satellites.append(Satellite(name, self.path, dataFile))
        
        
                
    