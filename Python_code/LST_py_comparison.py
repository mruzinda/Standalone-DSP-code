from datetime import datetime,timezone

# This is a test case to compare with calculate_phase function in C code
# (phase_calculation.c). It replicates the LST calculation in the Mosaic 
# python code.

year = 1969
month = 1
day = 5
ut = 25.083333333333333333333
inst_long = 20.4439

JD = 367*year - int(7*(year + int((month+9)/12))/4) - int(3*(int((year +\
 (month-9)/7)/100)+1)/4) + int(275*month/9) + day + 1721028.5 + ut/24 
print("Julian Date: ",JD) # For test -> 2440227.545139
J2000 = JD - 2451545.0
print("Days since J2000: ",J2000) # For test -> -11317.454861

# Local sidereal time (LST) as calculated in the Mosaic code
# with the following doc: http://www.stargazing.net/kepler/altaz.html#twig02
LST = 100.46 + 0.985647*J2000 + inst_long + 15*ut
print("LST before mod in degrees: ",LST) # For test -> -10657.861531

LST %= 360.0
print("LST in degrees: ",LST) # For test -> 142.138469

lst_hr = LST/360*24
lst_min = (lst_hr - int(lst_hr))*60
lst_sec = (lst_min - int(lst_min))*60
print(int(lst_hr), ":", int(lst_min), ":", int(lst_sec)) # For test -> 9:28:33