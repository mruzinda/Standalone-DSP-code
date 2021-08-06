from datetime import datetime,timezone

year = 1969
month = 1
day = 5

JD = 367*year - (7*(year + ((month+9)/12))/4) - (3*(((year + (month-9)/7)/100)+1)/4) +\
 (275*month/9) + day + 1721026.5
 
now_utc = datetime.now(timezone.utc)