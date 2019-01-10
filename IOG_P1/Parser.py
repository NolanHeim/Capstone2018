# -*- coding: utf-8 -*-
#
# Parser.py
#
# Parser has two main functions:
#
# 1. Parses "ephemeris" (.e) files with satellite data. This occurs when this file is run
#  standalone. Saves a numpy array in the 
#  Assumes the files have "EphemerisTimePosVel" before the satellite orbital data
#   and "END Ephemeris" at the end of the satellite data.
#  Assumes data is in the format:
#   Time (s), X (m), Y(m), Z(m), Vx(m), Vy(m), Vz(m)
#   Where X, Y, and Z are the positions and Vx, Vy and Vz are velocities.
#
# 2. Parses mission files (.txt) stored in the directory specified by missionpath. 
#  Creates corresponding mission objects that are returned as a list. This functionality
#  is not ran when the code is run standalone. The function to do this gets called within
#  the Mission Creator class.
# 
# Inputs:
#   datapath: The relative location of the satellite data. (ex: ../../Data/)
#   missionpath: The relative location of all satellite mission files (ex: ../../Missions/)
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
import os
import numpy as np
from Mission import *
import time
import uuid

###IMPORTANT - CHANGE THIS WHEN WE MOVE ON FROM PHASE 1###
phase1 = True
###


class Parser:
    
    def __init__(self, datapath):
        self.datapath = datapath        
        self.jDate_offset = 3
            
   
   # This function parses all required data from satellite .e files.
   # It outputs the orbital data (times, positions, velocities) and extra info
   # (JDate, UUID) to the disk, in the folder specified by parsed_datapath
    def parse_orbit_data(self, parsed_datapath):
        #for timing the parser
        t0 = time.time()
        
        for filename in os.listdir(self.datapath):
            # Open and parse the file
            dataFile = os.path.join(self.datapath, filename)
            fileObj = open(dataFile, 'r')
            fileLines = [line.rstrip('\n').strip(' ').strip('\t') for line in fileObj]
            fileLines[:] = [x for x in fileLines if x != '']
            fileObj.close()
        
            # Find the JDate in the file
            JDateIndex = fileLines.index('BEGIN Ephemeris') + self.jDate_offset
            JDLine = fileLines[JDateIndex].split()
            jDate = float(JDLine[-1])
        
            # Assign a UUID to this satellite.
            satelliteID = uuid.uuid4()        
 
            # Find the actual orbital data, generate the matrix
            dataStartIndex = fileLines.index('EphemerisTimePosVel') + 1
            dataStopIndex = fileLines.index('END Ephemeris') - 1
            dataMatrix = [[float(i) for i in fileLines[j].split()] for j in range(dataStartIndex,dataStopIndex)]
            
            # Path and filename for the orbital data matrix 
            newFile = os.path.join(parsed_datapath, 'parsed_'+filename[0:-2])

            # Path and filename for how the extra info will be stored.
            extraInfoFile = os.path.join(parsed_datapath, 'ex_parsed_'+filename[0:-2])
            
            
            # The zeros can later be replaced by other information as needed. 
            # For now these slots are unused. Only currently needed extra info
            # is the JDate (start of orbit data) and satellite ID
            extraInfo = [jDate, str(satelliteID), 0, 0, 0, 0, 0]            
            
            #save the parsed satellite data to the disk
            data = np.array(dataMatrix)
            np.save(newFile, data)
                        
            #save the extra info in the same way
            extra = np.array(extraInfo)
            np.save(extraInfoFile, extra)
            
        print("Total time = "+str(time.time() - t0))
        

# Code that runs the parser as standalone. This just parses and stores satellite data. 
if __name__ == "__main__":
    parser = Parser('../../Data/Orbit/')
    parser.parse_orbit_data('../../Parsed Data/')

