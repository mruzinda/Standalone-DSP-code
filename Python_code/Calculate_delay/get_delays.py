from katpoint import katpoint
import coordinate as coord
import os
import sys
import numpy as np
from coordinate import Antenna
import csv
import struct
from array import array
import time
from os.path import exists as file_exists

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
    def __init__(self, freq):
        """
        constructor of the Delay Polynomial class.
        
        """
        #A list of targets to observe
        targets = []

        # Beam idx 1
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        #print(target_string)
        targets.append(katpoint.Target(target_string))

        # Beam idx 2
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        #print(target_string)
        targets.append(katpoint.Target(target_string))

        # Beam idx 3
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        #print(target_string)
        targets.append(katpoint.Target(target_string))

        # Beam idx 4
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 5
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 6
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 7
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 8
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 9
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 10
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 11
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 12
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 13
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 14
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 15
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 16
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 17
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 18
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 19
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 20
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 21
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 22
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 23
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 24
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 25
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 26
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 27
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 28
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 29
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 30
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 31
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 32
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 33
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 34
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 35
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 36
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 37
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 38
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 39
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 40
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 41
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 42
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 43
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 44
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 45
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 46
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 47
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 48
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 49
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 50
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 51
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 52
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 53
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 54
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 55
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 56
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 57
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 58
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 59
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 60
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 61
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        #ra=":".join(['1','2','3'])
        #dec=":".join(['-3','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 62
        ra=":".join(['1','26','0'])
        dec=":".join(['-30','43','0'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 63
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        # Beam idx 64
        ra=":".join(['7','8','9'])
        dec=":".join(['-4','5','6'])
        target_string=",".join(['radec',ra,dec])
        targets.append(katpoint.Target(target_string))

        #Default boresight is at ("1:25:46.5336", "-30:42:39.815999")
        ra=":".join(['1','25','46.5336'])
        dec=":".join(['-30','42','39.815999'])
        target_string=",".join(['radec',ra,dec])
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
        self.frequency = freq #1.4e9
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
        #print(target_1d[64*2+2*2])

        return target_1d

    # Add method to return phase solutions from pickle file containing telstate phase solutions

def dict_to_antenna_ordered_list(dict_obj, antennas, pol='h'):
    ordered_list = []
    for antenna in antennas:
        antenna_key = "{}{}".format(antenna.name, pol)
        ordered_list.append(dict_obj[antenna_key])

    return ordered_list

# Compute delay polynomials and write to FIFO periodically
freq_tmp = 1.4e9
epoch_sec = 1629380016 #provide a random time for now
sim_flag = 0 # If sim_flag=1, then use simulated polynomials for testing
# Path to be created
path_c = "/tmp/obsfreq"
path_e = "/tmp/epoch"
path = "/tmp/katpoint_delays"
# If either one of these FIFOs exists, delete the files on startup
# os.unlink("pathtofile") -> Removes FIFO (remove file acting as FIFO)
if file_exists(path_c) == True:
    os.unlink(path_c)
if file_exists(path_e) == True:
    os.unlink(path_e)
if file_exists(path) == True:
    os.unlink(path)

blk_count = 0
obsfreq = 0

while 1:
    #--------Simulated data for debugging-------#
    if sim_flag == 1:
        output_tmp = np.zeros(8192)
        for i in range(0,8192):
            output_tmp[i] = float(i*1e-11)

        output1 = output_tmp.ravel().tolist()
        #print("Length of output1 array: ", len(output1))
        #print(type(output1[0]))
        #print(" ")
        #print(output1[0:15])
        #print(output1[128])
    #------------------------------------------#

    if blk_count == 0:
        blk_count = blk_count+1
        print(path_c)
        print(file_exists(path_c))
        while file_exists(path_c) == False and obsfreq == 0: # Assuming obsfreq will never be 0
            #print("In while file_exists == False")
            if file_exists(path_c) == True:
                print("In if file_exists == True")
                print(path_c)
                print(file_exists(path_c))
            
                with open(path_c, 'rb') as fifoc:
                    print("In with open(... rb) for fifoc")
                    obsfreq_tmp = struct.unpack('d', fifoc.read(8))
                    obsfreq = obsfreq_tmp[0]*1e6
                    print("Center frequency (Hz): ", obsfreq)
                os.unlink(path_c)
                break
            else:
                continue

        if file_exists(path_c) == True:
            print("In second if file_exists == True")
            print(path_c)
            print(file_exists(path_c))
            
            with open(path_c, 'rb') as fifoc:
                print("In with open(... rb) for fifoc")
                obsfreq_tmp = struct.unpack('d', fifoc.read(8))
                obsfreq = obsfreq_tmp[0]*1e6
                print("Center frequency (Hz): ", obsfreq)
            os.unlink(path_c)

        #print("Center frequency (Hz): ", obsfreq)

    print(path_e)
    print(file_exists(path_e))
    while file_exists(path_e) == False:
        #print("In while file_exists == False")
        if file_exists(path_e) == True:
            print("In if file_exists == True")
            print(path_e)
            print(file_exists(path_e))
            
            with open(path_e, 'rb') as fifo1:
                print("In with open(... rb) for fifo1")
                epoch_tmp = struct.unpack('d', fifo1.read(8))
                epoch_sec = epoch_tmp[0]
                print("Epoch in seconds: ", epoch_sec)
            os.unlink(path_e)
            break
        else:
            continue

    if file_exists(path_e) == True:
        print("In second if file_exists == True")
        print(path_e)
        print(file_exists(path_e))
            
        with open(path_e, 'rb') as fifo1:
            print("In with open(... rb) for fifo1")
            epoch_tmp = struct.unpack('d', fifo1.read(8))
            epoch_sec = epoch_tmp[0]
            print("Epoch in seconds: ", epoch_sec)
        os.unlink(path_e)

    print(path_e)
    print(file_exists(path_e))
    print("Epoch in seconds: ", epoch_sec)

    test = DelayPolynomial(obsfreq)
    output = test.get_delay_polynomials(epoch_sec,duration=2)

    print("Length of output array: ", len(output))
    print(type(output[0]))

    if file_exists(path) == True:
        os.unlink(path)
        print("Deleted katpoint_delays FIFO")

    try:
        os.mkfifo(path)
    except OSError as e:
        print("Failed to create FIFO: {0}" .format(e))
    else:
        with open(path, 'wb') as fifo:
            print("Before write to FIFO!")
            fifo.write(struct.pack('f'*len(output), *output))
            fifo.close()
            print("Path is created")

    # Remove FIFO (remove file acting as FIFO)
    os.unlink(path)

    # First beam set to boresight, so all delay should be zero
    # #print(output[0,:,0])
    # #Second beam 
    # #print(output[1,0:3,0])
    # #print(output[1,0:3,1])

    # The indices of the 1D arrays below match the ones commented above for the 3D array from get_delay_polynomials
    # First beam set to boresight, so all delay should be zero
    print(output[0:63])
    # Second beam
    print("[",output[64*2],output[64*2+2],output[64*2+2*2],"]")
    print("[",output[(64*2)+1],output[(64*2)+2+1],output[(64*2)+(2*2)+1],"]")

