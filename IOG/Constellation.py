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
    def __init__(self, datapath):
        self.path = datapath
        self.constellation = [];
        self.satellites = [];
        #initialize constellations
        
        
    #gets the info for one satellite
    def addToConstellation(self, satellite):
        self.constellation.append(satellite)

    def newSatellite(self, name, dataFile):
        self.satellites.append(Satellite(name, self.path, dataFile))
        
        
                
    