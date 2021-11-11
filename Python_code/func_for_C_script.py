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
    def __init__(self):
        """
        constructor of the Vals_to_generate class.

        """
        self.array_vals = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]

    def gen_some_vals(self):
        """
        Return the values of the array initialized in __init__()

        """
        return self.array_vals

#p1 = Vals_to_generate()
#x = p1.gen_some_vals()

#print(x[3])
#print(x)

