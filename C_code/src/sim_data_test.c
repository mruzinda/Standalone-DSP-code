#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "coherent_beamformer_char_in.h"

// Used to investigate the simulated data a little further.

int main(){
	signed char* data_sim;
	data_sim = (signed char*)calloc(N_INPUT, sizeof(signed char));

	float c = 3e8; // Speed of light
	float c_freq = 1.25e9; // Center frequency
	float lambda = c / c_freq; // Wavelength
	float d = lambda / 2; // Distance between antennas
	float chan_band = 1; // Fine channel bandwidth in Hz

	//float* rf_freqs = (float*)calloc(N_FREQ, sizeof(float));
	//for (int i = 0; i < N_FREQ; i++) {
	//	rf_freqs[i] = chan_band * i + c_freq;
	//}

	//float* theta = (float*)calloc(N_TIME, sizeof(float)); // SOI direction/angle of arrival
	//float* tau = (float*)calloc(N_TIME, sizeof(float)); // Delay

	float theta = 0; // SOI direction/angle of arrival
	float tau = 0; // Delay
	float rf_freqs = 0;
	float tmp_max = 1.0;
	float tmp_min = -1.0;

	for (int t = 0; t < N_TIME; t++) {
		theta = ((t/100 - (N_TIME / 200)) + 90)*PI/180; // SOI direction/angle of arrival -> Moving across array over time i.e. angle changes each time sample
		tau = d * cos(theta) / c; // Delay
		for (int f = 0; f < N_FREQ; f++) {
			rf_freqs = chan_band * f + c_freq;
			for (int a = 0; a < N_ANT; a++) {
				if(a < N_REAL_ANT){
					// X polarization
					data_sim[2 * data_in_idx(a, 0, f, t)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min))-0.5)*256);
					//data_sim[2 * data_in_idx(a, 0, f, t)] = (signed char)(((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
					data_sim[2 * data_in_idx(a, 0, f, t) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min))-0.5)*256);
					// Y polarization
					data_sim[2 * data_in_idx(a, 1, f, t)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min))-0.5)*256);
					data_sim[2 * data_in_idx(a, 1, f, t) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min))-0.5)*256); // Make this negative if a different polarization is tested
					if((a >= 0) && (a < 6) && (f == 2) && (t == 0)){
						printf("X - data_sim r: %d, data_sim i: %d\n", data_sim[2 * data_in_idx(a, 0, f, t)], data_sim[2 * data_in_idx(a, 0, f, t) + 1]);
						printf("Y - data_sim r: %d, data_sim i: %d\n", data_sim[2 * data_in_idx(a, 1, f, t)], data_sim[2 * data_in_idx(a, 1, f, t) + 1]);
						printf("Y - data_in_idx: %d\n", 2 * data_in_idx(a, 0, f, t));
						printf("cos(): %lf\n", cos(2 * PI * rf_freqs * a * tau));
						printf("sin(): %lf\n", sin(2 * PI * rf_freqs * a * tau));
					}
				}else{
					// X polarization
					data_sim[2 * data_in_idx(a, 0, f, t)] = 0;
					data_sim[2 * data_in_idx(a, 0, f, t) + 1] = 0;
					// Y polarization
					data_sim[2 * data_in_idx(a, 1, f, t)] = 0;
					data_sim[2 * data_in_idx(a, 1, f, t) + 1] = 0; // Make this negative if a different polarization is tested
				}
			}
		}
	}
	for (int a = 0; a < N_ANT; a++){
		printf("Antenna %d = %d\n", (a+1), data_sim[2 * data_in_idx(a, 0, 0, 0)]);
	}
	return 0;
}
