from datetime import datetime,timezone

year = 1969
month = 1
day = 5

JD = 367*year - (7*(year + ((month+9)/12))/4) - (3*(((year + (month-9)/7)/100)+1)/4) +\
 (275*month/9) + day + 1721026.5
print("Julian Date: ",JD)
J2000 = JD - 2451545.0
print("Days since J2000: ",J2000)

ut = 25.083333333333333333333
inst_long = 20.4439
LST = 100.46 + 0.985647*J2000 + inst_long + cur_time/24
print("LST before mod in degrees: ",LST)

LST %= 360.0
print("LST in degrees: ",LST)

lst_hr = LST/360*24
lst_min = (lst_hr - int(lst_hr))*60
lst_sec = (lst_min - int(lst_min))*60
print(int(lst_hr), ":", int(lst_min), ":", int(lst_sec))