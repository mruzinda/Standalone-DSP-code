from katpoint import katpoint
import coordinate as coord
import sys
import numpy as np
from coordinate import Antenna
import csv

#sys.path.append('/home/mruzinda/Calculate_delay/katpoint/katpoint')

#import delay

class DelayPolynomial(object):
    """
    Class for generation of  delay polynomial
    
    arguments:
    antennas -- a list of antenna objects or coordinate in csv format
    targets -- a list of beam location in equatorial coordinates
    frequencies -- a list of frequencies on which the polynomail is calculated in Hz
    reference -- the reference antenna for delay calculation
    """
    #def __init__(self, antennas, bore_sight, targets, reference):
    def __init__(self):
        """
        constructor of the Delay Polynomial class.
        
        """
        #A list of targets to observe
        targets = []

        #ha=":".join(['1','25','46.5336'])
        #dec=":".join(['-30','42','39.815999'])
        ha=":".join(['1','2','3'])
        dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ha,dec])
        print(target_string)
        targets.append(katpoint.Target(target_string))

        ha=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ha,dec])
        print(target_string)
        targets.append(katpoint.Target(target_string))

        ha=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ha,dec])
        print(target_string)
        targets.append(katpoint.Target(target_string))

        #Default boresight is at ("1:25:46.5336", "-30:42:39.815999")
        ha=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        target_string=",".join(['radec',ha,dec])
        print(target_string)
        #defaultBoreSight=katpoint.Target(target_string)
        bore_sight=katpoint.Target(target_string)

        #Reference antenna
        arrayRef = katpoint.Antenna('ref, -30:42:39.8, 21:26:38.0, 1035.0')
        arrayRef.position_wgs84

        #Load in antennas 
        """
        Need a list of antenna objects or coordinate in csv format
        katpoint.Antenna class says that the format for a fully-specified antenna should be: 
        'FF2, -30:43:17.3, 21:24:38.5, 1038.0, 12.0, 86.2 25.5 0.0, -0:06:39.6 0 0 0 0 0 0:09:48.9, 1.16'
        
        """
        ants=[]
        with open('/home/mruzinda/Calculate_delay/antenna.csv', newline='') as csvfile:
            data = csv.reader(csvfile, delimiter=',')
            for row in data:
                ants.append(katpoint.Antenna(', '.join(row)))

        #self.antennas = antennas
        self.antennas = ants
        self.targets = DelayPolynomial.check_targets(targets)
        self.frequency = 1.4e9
        #self.reference = reference
        self.reference = arrayRef
        self.bore_sight = DelayPolynomial.check_targets([bore_sight,])[0]

    @staticmethod
    def check_targets(targets):
        """
        check the target data type, the arguments will be converted to katpoint
            object if they are not.
        
        arguments:
        targets -- a list of target objets in the format of
            katpoint target object or set of [longitude, latitude, altitude]
        
        return:
            targets in katpoint object
        
        """
        if isinstance(targets[0], katpoint.Target):
            return targets
        else:
            return DelayPolynomial.make_katpoint_target(targets)

    @staticmethod
    def check_time(time):
        """
        check the the data type of the time value. If the values are datetime
            objects, they will be converted to seconds
        
        arguments:
        time -- epoch seconds or datetime objects
        
        return:
        time in epoch seconds
        """
        if type(time) != int and type(time) != float:
            return coord.datetimeToEpoch(time)
        else:
            return time

    @staticmethod
    def make_katpoint_target(sources):
        """
        check the the data type of the source. If the values are in (RA, DEC) pair,
        they will be converted to katpoint target object
        
        arguments:
        source -- source in either (RA, DEC) pairs or katpoint target objects
        
        return:
        katpoint target objects
        """

        targets = []
        for source in sources:
            target_string = ",".join(['radec',
                            coord.angleToHour(source[0]),
                            coord.angleToDEC(source[1])])
            targets.append(katpoint.Target(target_string))
        return targets

    def get_delay_polynomials(self, epoch, duration=10.0):
        """
        calculate and return the polynomials
        
        Arguments:
        timestamp -- the observation time in datatime object or epoch seconds
        duration -- the duration in which the polynomial is calcuated
        
        return:
        polynomials in the order of beam, antenna, (delay, rate)
        
        """
        timestamp = DelayPolynomial.check_time(epoch)
        antennaObjectList = self.antennas
        timestamp = (timestamp, timestamp + duration)

        target_array = []
        for target in self.targets:
            dc = katpoint.DelayCorrection(self.antennas,
                    self.reference, self.frequency)
            delay, phase = dc.corrections(target, timestamp)
            delayArray = np.array(dict_to_antenna_ordered_list(
                        delay, self.antennas))

            target_array.append(delayArray[:,0,:])
        target_array = np.array(target_array)

        dc = katpoint.DelayCorrection(self.antennas,
            self.reference, self.frequency)
        delay, phase = dc.corrections(self.bore_sight, timestamp)
        bore_sight_delay = np.array(dict_to_antenna_ordered_list(
                    delay, self.antennas))[:,0,:]

        target_array = target_array - bore_sight_delay
        # return target_array
        dim_1 = len(target_array[:,0,0])
        dim_2 = len(target_array[0,:,0])
        dim_3 = len(target_array[0,0,:])
        target_1d_tmp = target_array.reshape(1,dim_1*dim_2*dim_3)
        target_1d = target_1d_tmp.ravel().tolist()
        print(target_1d[64*2+2*2])

        return target_1d

def dict_to_antenna_ordered_list(dict_obj, antennas, pol='h'):
    ordered_list = []
    for antenna in antennas:
        antenna_key = "{}{}".format(antenna.name, pol)
        ordered_list.append(dict_obj[antenna_key])

    return ordered_list

# #provide a random time for now
# test = DelayPolynomial()
# time = 1629380016
# output = test.get_delay_polynomials(time,duration=2)

# #beam, antenna, (delay, rate)
# print(output.shape)

# #First beam set to boresight, so all delay should be zero
# #print(output[0,:,0])
# print(output[0,0:63])

# #Second beam 
# #print(output[1,0:3,0])
# #print(output[1,0:3,1])
# # The indices of the 1D arrays below match the ones commented above
# print(output[0,[64*2,64*2+2,64*2+2*2]])
# print(output[0,[(64*2)+1,(64*2)+2+1,(64*2)+(2*2)+1]])




