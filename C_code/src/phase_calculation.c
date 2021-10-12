#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/types.h>
//#include <fftw3.h>
#include <unistd.h>
//#include "hashpipe.h"
//#include "meerkat_databuf.h"
#include "metadata.h" 

#define PI 3.14159265
#define light 299792458.
#define one_over_c 0.0033356
#define R2D (180. / PI)
#define D2R (PI / 180.)
#define TAU (2 * PI)
#define inst_long  20.4439 // 21.23
#define inst_lat -30.7111  // 30.42  //Do we need to take altitude into account?
#define test_lst 0 // Used to check whether LST is being calculated correctly
#define VERBOSE 1 // Print values to terminal from function

void calculate_phase(int nbeams, int nants, uint32_t nchans_in, struct timespec time_now, float* freq_now, struct beamCoord beam_coord, struct antCoord ant_coord, double* phase)
{
  //struct tm* timeinfo;
  //timeinfo = localtime(&time_now.tv_sec);
  time_t timeNow = time(0);
  struct tm* timeinfo = localtime(&timeNow);
  int year = timeinfo->tm_year + 1900;
  int month = timeinfo->tm_mon + 1;
  /*
  if (month < 3) {
    month = month + 12;
    year = year - 1;
  }
  */
  int day = timeinfo->tm_mday;
  
  #if VERBOSE == 1
  printf("Time info: year = %d, month = %d, day = %d, hour = %d, minutes = %d\n", timeinfo->tm_year, timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min);
  printf("Year = : %d, month = %d, day = %d\n", year, month, day);
  #endif
  
  //float JD = 2 - (int)(year / 100.) + (int)(int)((year / 100.) / 4.) + (int)(365.25 * year)
  //  + (int)(30.6001 * (month + 1)) + day + 1720994.5 + 0.5; //+ 1720994.5; // JD - whole days since last epoch
  //printf("JD: %lf\n", JD);
  //double T = (JD - 2451545.0)
  //  / 36525.0; // Works if time after year 2000, otherwise T is -ve and might break
  //printf("T: %lf\n", T);
  ////double T0 = fmod((6.697374558 + (2400.051336 * T) + (0.000025862 * T * T)), 24.);
  //double T0 = 6.697374558 + (0.06570982441908*JD) + (0.000025862 * T * T);
  //double UT = (timeinfo->tm_hour) + (timeinfo->tm_min / 60.)
  //  + (timeinfo->tm_sec + time_now.tv_nsec / 1.e9) / 3600.; // Hours since midnight - Correct
  //printf("UT: %lf\n", UT); // Correct
  //double GST = fmod((T0 + UT * 1.002737909), 24.);
  //double GST = 24110.54841 + 8640184.812866*T + 0.093104*T*T - 0.0000062*T*T*T;
 
  //printf("GST: %lf\n", GST);
  //double LST = GST + inst_long; // / 15.;
  //while (LST < 0) {
  //  LST = LST + 24;
  //}
  //LST = fmod(LST, 24);
  
  double UT = (timeinfo->tm_hour) + (timeinfo->tm_min / 60.)
    + (timeinfo->tm_sec + time_now.tv_nsec / 1.e9) / 3600.; // Hours since midnight - Correct
  
  #if VERBOSE == 1
  printf("UT: %lf\n", UT); // Correct
  #endif

  // This is a test case to compare with python code used by Mosaic 
  // (another project on MeerKAT). It has fixed values that help compare 
  // the LST with the LST_py_comparison.py script which replicates the 
  // Mosaic LST calculation in python.
  #if test_lst == 1
  year = 1969;
  month = 1;
  day = 5;
  UT = 25.083333333333333;
  #endif
  
  // Julian Date calculation
  double JD = 367*year - (7*(year + ((month+9)/12))/4) - (3*(((year + (month - 9)/7)/100)+1)/4) 
	+ (275*month/9) + day + 1721028.5 + UT/24;
  
  #if VERBOSE == 1
  printf("JD: %lf\n", JD); // For test -> 2440227.545139
  #endif
  
  // Days since J2000
  double T = (JD - 2451545.0); 
  
  #if VERBOSE == 1
  printf("T: %lf\n", T); // For test -> -11317.454861
  #endif
  
  // Local sidereal time (LST) as calculated in the Mosaic code
  // with the following doc: http://www.stargazing.net/kepler/altaz.html#twig02
  double LST = 100.46 + 0.985647*T + inst_long + 15*UT;
  
  #if VERBOSE == 1
  printf("LST before mod by 360: %lf\n", LST); // For test -> -10657.861531
  #endif

  // LST in degrees used for the rest of the calculation
  // rather than in hours
  LST = remainder(LST,360);

  #if VERBOSE == 1
  printf("LST in degrees: %lf\n", LST); // For test -> 142.138469
  #endif
  
  #if test_lst == 1
  // LST in hours, minutes, and seconds
  float LST_hr = LST/360*24; 
  float LST_min = (LST_hr - (int)LST_hr)*60;
  float LST_sec = (LST_min - (int)LST_min)*60;

  // For test -> 9:28:33
  printf("LST in hours:min:sec = %d:%d:%d\n", (int)LST_hr, (int)LST_min, (int)LST_sec);
  printf("LST (float) in hours:min:sec = %f:%f:%f\n", LST_hr, LST_min, LST_sec);
  #endif
  // LST in hours, minutes, and seconds
  float LST_hr = LST/360*24; 
  float LST_min = (LST_hr - (int)LST_hr)*60;
  float LST_sec = (LST_min - (int)LST_min)*60;

  printf("LST in hours:min:sec = %d:%d:%d\n", (int)LST_hr, (int)LST_min, (int)LST_sec);

  for (int b = 0; b < nbeams; b++) {
    double hour_angle = LST * 15. - beam_coord.ra[b];
    double alt = sin(beam_coord.dec[b] * D2R) * sin(inst_lat * D2R)
                 + cos(beam_coord.dec[b] * D2R) * cos(inst_lat * D2R) * cos(hour_angle * D2R);
    alt = asin(alt);
    double az = (sin(beam_coord.dec[b] * D2R) - sin(alt) * sin(inst_lat * D2R))
      / (cos(alt) * cos(inst_lat * D2R));
    az = acos(az);
    if (sin(hour_angle * D2R) >= 0) {
      az = TAU - az;
    }
    double projection_angle, effective_angle, offset_distance;
    for (int i = 0; i < nants ; i++){
      float dist_y = ant_coord.north[i];
      float dist_x = ant_coord.east[i];
      projection_angle = 90 * D2R - atan2(dist_y, dist_x); //Need to check this
      offset_distance = sqrt(pow(dist_y, 2) + pow(dist_x, 2));
      effective_angle = projection_angle - az;
      for (int f = 0; f < nchans_in ; f++ ) {
	//phase [freq-beam-ant]
	int ph_id = f*nbeams*nants + b*nants + i;
	phase[ph_id*2] = cos(TAU * cos(effective_angle) * cos(alt) * offset_distance
			     * freq_now[f] * one_over_c);
	phase[ph_id*2+1] = -sin(TAU * cos(effective_angle) * cos(-alt) * offset_distance
				* freq_now[f] * one_over_c);
      }
    }
  }
}

int main(){

    int nbeams=64;
    int nants=64;
    int npols=2;
    int nupchan = 4; //should be 2**19 for 1k mode
    int nchans_in = 1;
    int nchans_out = nchans_in*nupchan; 
    int nsamps_in = 8*nupchan;
    int nsamps_out = nsamps_in/nupchan;

    int input_len = nants*nchans_in*nsamps_in*npols*2;
    int output_len = nbeams*nchans_out*nsamps_out;
    float* host_output;
    char* host_input;
    double* host_phase;
    float* upchan_output;
    host_output = (float*)malloc(output_len * sizeof(float));
    host_input = (char*)malloc(input_len * sizeof(char));
    host_phase = (double*)malloc(nchans_in*(nbeams-1)*nants*2*sizeof(double)); //assume one phase correction for each coarse freq
    upchan_output = (float*)malloc(input_len*sizeof(float));
    float* freq_now = (float*)malloc(nchans_in*sizeof(float));
    struct timespec time_now;
    struct beamCoord beam_coord;
    struct antCoord ant_coord;
    bool update_antenna = true;

    //Freq
    for (int f = 0 ; f < nchans_in ; f++ ) {
        freq_now[f] = 700.0 + f*0.00836;
    }

    for (int b = 0; b< nbeams-1; b++){
	beam_coord.ra[b] = 128.83588121;
        beam_coord.dec[b] = -45.17635419;
    }

    //Antenna position (read from lookup table, only need to do it once)----------
    if (update_antenna) {
        FILE* ptr = fopen("meerkat_antenna_positions.dat","r");
        if (ptr==NULL)
        {
            printf("Missing antenna positions");
	    return 0;
        }
        fscanf(ptr, "%*[^\n]\n"); //skip first line
        for (int i = 0 ; i < nants; i++ ) {
            if (fscanf(ptr, "%*s %f %f %f", &ant_coord.north[i], &ant_coord.east[i], &ant_coord.up[i]) != 3){
	        printf("problem reading in positions");
	        return 0;
	    }
        }
        fclose(ptr);
        update_antenna  = false;
    }

    //Calcuate phase------------------------------------------------------------
    calculate_phase(nbeams-1, nants, nchans_in, time_now, freq_now, beam_coord, ant_coord, host_phase);

    printf("Array response of beam 1 element 1; Real = %lf and Imag = %lf", host_phase[0], host_phase[1]);

    return 0;
}