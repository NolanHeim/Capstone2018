# -*- coding: utf-8 -*-
#
# CLASSNAME.py
#
# DESCRIPTION
#
# Initial Creation Date: 09/26/2018
#
# Written by Jordan Jones and Nolan Heim
#
#For numpy.piecewise
import numpy

class Calculator:
    
    def __init__(self):
        self.constant = 5
    
    #Returns the cubic Hermite polynomial function on the
    #subinterval [subt1,subt2]
    def cubicHermitePolynomial(self, subt1, subt2, func):
        condition = lambda t: subt1 <= t <= subt2        
        h = subt2-subt1
        f1 = func[subt1]
        f2 = func[subt2]
        Ci = lambda t: h*f1*f2*t#Function Definition Here
        
        return [condition, Ci]
        
    def visibilityFunction(self, t1, t2, func, N):
        print('Xi')