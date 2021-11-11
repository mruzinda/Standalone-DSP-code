# Python script with functions for test C script to call

import matplotlib.pyplot as plt
from array import array
import numpy as np
import sys

class Vals_to_generate(object):
    """
    Generates values for testing module call from C script.
    Or it can be used for other purposes as well.
    """
    def __init__(self,func_flag):
        """
        Constructor of the Vals_to_generate class.
        Argument is just a flag to alter the data and test argument input

        """
        self.func_flag = func_flag
	if self.func_flag == 0:
            #self.array_vals = 2
            self.array_vals = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        elif self.func_flag == 1:
            #self.array_vals = 3
            self.array_vals = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

    def gen_some_vals(self):
        """
        Return the values of the array initialized in __init__()

        """
        return self.array_vals

#p1 = Vals_to_generate(1)
#x = p1.gen_some_vals()

#print(x[0])
#print(x)

