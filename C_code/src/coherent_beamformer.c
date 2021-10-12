#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "coherent_beamformer.h"

// Generate simulated data
float* simulate_data();

// Generate simulated weights or coefficients
float* simulate_coefficients();

// Perform transpose on the data and convert to floats
void data_transpose(float* data_in, float* data_tra);

// Perform beamforming operation
void coherent_beamformer(float* input_data, float* coeff, float* output_data);

// Compute power of beamformer output (abs()^2)
void beamformer_power(float* bf_volt, float* bf_power);

//float* h_data = NULL;
float* data_float = NULL;
float* coh_bf_out = NULL;
// Allocate memory to all arrays 
void init_beamformer() {
	// Allocate memory for input data float type for transpose
	data_float = (float*)calloc(N_INPUT, sizeof(float));
	// Allocate memory for coherent beamformer output
	coh_bf_out = (float*)calloc(N_OUTPUT, sizeof(float));

	return;
}

// Perform transpose on the data and convert to floats
void data_transpose(float* data_in, float* data_tra) {
	// a - Antenna index
	// p - Polarization index
	// f - Frequency index
	// t - Time sample index

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	for (int t = 0; t < N_TIME; t++){
		for (int f = 0; f < N_FREQ; f++){
			for (int p = 0; p < N_POL; p++){
				for (int a = 0; a < N_ANT; a++){
					data_tra[2*data_tr_idx(a, p, f, t)] = data_in[2*data_in_idx(a, p, f, t)];
					data_tra[2*data_tr_idx(a, p, f, t)+1] = data_in[2*data_in_idx(a, p, f, t) + 1];
				}
			}
		}
	}
	
	return;
}

// Perform beamforming operation
void coherent_beamformer(float* input_data, float* coeff, float* output_data) {
	// a - Antenna index
	// p - Polarization index
	// f - Frequency index
	// t - Time sample index
	// b - Beam index

	float bf_product_re;
	float bf_product_im;
	float bf_in_data_re;
	float bf_in_data_im;
	float bf_coeff_re;
	float bf_coeff_im;

	for (int t = 0; t < N_TIME; t++){
		for (int f = 0; f < N_FREQ; f++){
			for (int b = 0; b < N_BEAM; b++){
				for (int p = 0; p < N_POL; p++){
					for (int a = 0; a < N_ANT; a++){
						int i = data_tr_idx(a, p, f, t);
						int w = coeff_idx(a, p, b, f);
						int h = coh_bf_idx(p, b, f, t);
						
						bf_in_data_re = input_data[2*i];
						bf_in_data_im = input_data[2*i + 1];
						bf_coeff_re = coeff[2*w];
						bf_coeff_im = coeff[2*w + 1];
						
						bf_product_re = (bf_in_data_re * bf_coeff_re) + (bf_in_data_im * bf_coeff_im);
						bf_product_im = (bf_in_data_re * bf_coeff_im) - (bf_in_data_im * bf_coeff_re);
						
						// Beamform (Sum all antennas)
						output_data[2*h] += bf_product_re;
						output_data[2*h + 1] += bf_product_im;
					}
				}
			}
		}
	}
	return;
}

// Compute power of beamformer output (abs()^2)
void beamformer_power(float* bf_volt, float* bf_power) {
	// p - Polarizaton index
	// b - Beam index
	// f - Frequency bin index
	// t - Time sample index

	// Power = Absolute value squared of output -> r^2 + i^2
	// And sum polarizations
	int xp;
	int yp;
	float x_pol_pow;
	float y_pol_pow;
	for (int t = 0; t < N_TIME; t++){
		for (int f = 0; f < N_FREQ; f++){
			for (int b = 0; b < N_BEAM; b++){
				xp = coh_bf_idx(0, b, f, t); // X polarization
				yp = coh_bf_idx(1, b, f, t); // Y polarization
				
				x_pol_pow = (bf_volt[2*xp]*bf_volt[2*xp]) + (bf_volt[2*xp + 1]*bf_volt[2*xp + 1]); // XX*
				y_pol_pow = (bf_volt[2*yp]*bf_volt[2*yp]) + (bf_volt[2*yp + 1]*bf_volt[2*yp + 1]); // YY*
				
				bf_power[pow_bf_idx(b, f, t)] = x_pol_pow + y_pol_pow; // XX* + YY*
			}
		}
	}
	
	
	return;
}

// Run beamformer
void run_beamformer(float* data_in, float* bf_coefficient, float* data_out) {

	float* data_tra = data_float;
	float* bf_output = coh_bf_out;

	// Perform transpose on the data and convert to floats  
	data_transpose(data_in, data_tra);
	
	// Perform beamforming operation
	coherent_beamformer(data_tra, bf_coefficient, bf_output);

	// Compute power of beamformer output (abs()^2)    
	beamformer_power(bf_output, data_out);

	return;
}

// Generate simulated data
float* simulate_data() {
	float* data_sim;
	data_sim = (float*)calloc(N_INPUT, sizeof(float));
	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Repeating sequence of 1 to 64
	sim_flag = 2 -> Sequence of 1 to 64 placed in a particular bin (bin 3 and 6 for now)
	sim flag = 3 -> Simulated radio source in center beam assuming ULA
	*/
	int sim_flag = 3;
	if (sim_flag == 0) {
		for (int i = 0; i < (N_INPUT / 2); i++) {
			data_sim[2 * i] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int t = 0; t < N_TIME; t++) {
				for (int f = 0; f < N_FREQ; f++) {
					for (int a = 0; a < N_ANT; a++) {
						if (tmp >= N_ANT) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_ANT+1);
						data_sim[2 * data_in_idx(a, p, f, t)] = tmp;
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int t = 0; t < N_TIME; t++) {
				for (int a = 0; a < N_ANT; a++) {
					if (tmp >= N_ANT) {
						tmp = 0;
					}
					tmp = (tmp + 1) % (N_ANT+1);
					data_sim[2 * data_in_idx(a, p, 5, t)] = tmp;
					data_sim[2 * data_in_idx(a, p, 2, t)] = tmp;
				}
			}
		}
	}
	if (sim_flag == 3) {
		float c  = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c/c_freq; // Wavelength
		float d = lambda/2; // Distance between antennas
		float chan_band = 1; // Fine channel bandwidth in Hz
		
		float * rf_freqs = (float *)calloc(N_FREQ, sizeof(float));
		for(int i = 0; i < N_FREQ; i++){
			rf_freqs[i] = chan_band*i + c_freq;
		}
		
		float * theta = (float *)calloc(N_TIME, sizeof(float)); // SOI direction/angle of arrival
		float * tau = (float *)calloc(N_TIME, sizeof(float)); // Delay
		
		for (int t = 0; t < N_TIME; t++) {
			theta[t] = (t - (N_TIME/2)) + 90; // SOI direction/angle of arrival -> Moving across array over time i.e. angle changes each time sample
			tau[t] = d*cos(theta[t])/c; // Delay
			for (int f = 0; f < N_FREQ; f++) {
				for (int a = 0; a < N_ANT; a++) {
					// X polarization
					data_sim[2 * data_in_idx(a, 0, f, t)] = cos(2 * PI * rf_freqs[f] * a * tau[t]);
					data_sim[2 * data_in_idx(a, 0, f, t) + 1] = sin(2 * PI * rf_freqs[f] * a * tau[t]);
					// Y polarization
					data_sim[2 * data_in_idx(a, 1, f, t)] = cos(2 * PI * rf_freqs[f] * a * tau[t]);
					data_sim[2 * data_in_idx(a, 1, f, t) + 1] = sin(2 * PI * rf_freqs[f] * a * tau[t]); // Make this negative if a different polarization is tested
				}
			}
		}
	}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients() {
	float* coeff_sim;

	coeff_sim = (float*)calloc(N_COEFF, sizeof(float));
	
	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64
	sim_flag = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2
	sim flag = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.
	*/
	int sim_flag = 3;
	if (sim_flag == 0) {
		for (int i = 0; i < (N_COEFF/2); i++) {
			coeff_sim[2*i] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int a = 0; a < N_ANT; a++) {
				for (int f = 0; f < N_FREQ; f++) {
					for (int b = 0; b < N_BEAM; b++) {
						if (tmp >= N_BEAM) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_BEAM+1);
						coeff_sim[2 * coeff_idx(a, p, b, f)] = tmp;
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int a = 0; a < N_ANT; a++) {
				for (int b = 0; b < N_BEAM; b++) {
					if (tmp >= N_BEAM) {
						tmp = 0;
					}
					tmp = (tmp + 1) % (N_BEAM+1);
					coeff_sim[2 * coeff_idx(a, p, b, 2)] = tmp;
					coeff_sim[2 * coeff_idx(a, p, b, 5)] = tmp;
				}
			}
		}
	}
	if (sim_flag == 3) {
		float c  = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c/c_freq; // Wavelength
		float d = lambda/2; // Distance between antennas
		float chan_band = 1; // Fine channel bandwidth in Hz
		
		float * rf_freqs = (float *)calloc(N_FREQ, sizeof(float));
		for(int i = 0; i < N_FREQ; i++){
			rf_freqs[i] = chan_band*i + c_freq;
		}
		
		float * theta = (float *)calloc(N_TIME, sizeof(float)); // Beam angle from 58 to 122 degrees
		float * tau_beam = (float *)calloc(N_BEAM, sizeof(float)); // Delay
		
		for (int b = 0; b < N_BEAM; b++) {
			theta[b] = (b - (N_BEAM/2)) + 90; // Beam angle from 58 to 122 degrees - Given SOI at 90 deg or moving across array, the beam with the most power is beamm 33
			tau_beam[b] = d*cos(theta[b])/c; // Delay
			for (int f = 0; f < N_FREQ; f++) {
				for (int a = 0; a < N_ANT; a++) {
					// X polarization
					coeff_sim[2 * coeff_idx(a, 0, b, f)] = cos(2 * PI * rf_freqs[f] * a * tau_beam[b]);
					coeff_sim[2 * coeff_idx(a, 0, b, f) + 1] = sin(2 * PI * rf_freqs[f] * a * tau_beam[b]);
					// Y polarization
					coeff_sim[2 * coeff_idx(a, 1, b, f)] = cos(2 * PI * rf_freqs[f] * a * tau_beam[b]);
					coeff_sim[2 * coeff_idx(a, 1, b, f) + 1] = sin(2 * PI * rf_freqs[f] * a * tau_beam[b]); // Make this negative if a different polarization is tested
				}
			}
		}
	}
	
	/* // Check to see whether these are the values expected
	for (int b = 0; b < N_BEAM; b++){
		printf("Beam %d = %f\n", (b+1), coeff_sim[2 * coeff_idx(20, 0, b, 2)]);
	}
	*/

	return coeff_sim;
}

// Free memory
void cohbfCleanup() {
	// Free up memory at the end of a program
	if (data_float != NULL) {
		free(data_float);
	}
	if (coh_bf_out != NULL) {
		free(coh_bf_out);
	}
}

// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");
	// Generate simulated data
	float* sim_data = simulate_data();

	printf("Here1!\n");

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients();

	printf("Here2!\n");

	// Allocate memory to all arrays used by run_beamformer() 
	init_beamformer();

	printf("Here3!\n");

	// Allocate memory for output array
	float* output_data;
	output_data = (float*)calloc(N_BF_POW, sizeof(float));

	printf("Here4!\n");

	// Run beamformer 
	run_beamformer(sim_data, sim_coefficients, output_data);

	printf("Here5!\n");

	int plt_flag = 0; // To plot figure here with gnuplot, set flag to 1
	
	if(plt_flag == 1){
		// Commands used for gnuplot
		char * commandsForGnuplot[] = {"set title \"Beamformer output\"", "plot 'intensitymap.txt' matrix with image"};
	
		// Plot intensity map of a single beam of the output power using gnuplot
		FILE *fp=NULL;
		fp=fopen("intensitymap.txt","w");
	
		float t0,t1,t2,t3,t4,t5,t6,t7; // Time samples 
		int b = 0; // Beam index
		for(int f=0; f<N_FREQ; f++){
			t0 = output_data[pow_bf_idx(b, f, 0)];
			t1 = output_data[pow_bf_idx(b, f, 1)];
			t2 = output_data[pow_bf_idx(b, f, 2)];
			t3 = output_data[pow_bf_idx(b, f, 3)];
			t4 = output_data[pow_bf_idx(b, f, 4)];
			t5 = output_data[pow_bf_idx(b, f, 5)];
			t6 = output_data[pow_bf_idx(b, f, 6)];
			t7 = output_data[pow_bf_idx(b, f, 7)];
			fprintf(fp,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", t0,t1,t2,t3,t4,t5,t6,t7);
		}
	
		FILE * gnuplotPIpe = popen ("gnuplot -persistent", "w");
  
		//fprintf(gnuplotPIpe, "plot 'heatmap.txt' matrix with image \n");
	  
		for(int i=0; i<N_COMMANDS; i++){
			fprintf(gnuplotPIpe, "%s \n", commandsForGnuplot[i]);
		}
	
		fflush(gnuplotPIpe);
	}
	
	// Write data to text file for analysis
	char output_filename[128];

	printf("Here6!\n");

	//strcpy(output_filename, "C:\Users\ruzie\OneDrive\Desktop\Work\CUDA_code\output_d.txt");
	strcpy(output_filename, "output_d_c.txt");

	printf("Here7!\n");

	FILE* output_file;

	printf("Here8!\n");

	output_file = fopen(output_filename, "w");

	printf("Here9!\n");

	for (int ii = 0; ii < N_BF_POW; ii++) {
		fprintf(output_file, "%g\n", output_data[ii]);
	}

	printf("Here10!\n");

	fclose(output_file);

	printf("Closed output file.\n");

	free(output_data);

	printf("Freed output array memory.\n");

	// Free up device memory
	cohbfCleanup();

	printf("Here11!\n");

	return 0;
}