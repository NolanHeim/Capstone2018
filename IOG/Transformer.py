# -*- coding: utf-8 -*-
#
# Transformer.py
#
# Initial Creation Date: 09/30/2018
#
# Written by Jordan Jones and Nolan Heim
#

import math

class Transformer:
    
    def ___init__(self):
        self.r_Earth = 6378137.0 #m
   
    
    def geo_2_eci(self, lat, lon, alt, theta): #theta in radians
        phi = theta + math.radians(lon)
        
        radius = alt + self.r_Earth*math.cos(math.radians(lat))
        x = radius*math.cos(phi)
        y = radius*math.sin(phi)
        z = (alt + self.r_Earth)*math.sin(math.radians(lat))

