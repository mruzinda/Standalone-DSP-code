#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <curand.h>
#include <assert.h>
//#include <unistd.h>
#include <cublas_v2.h>
#include <time.h>
//#include <sys/time.h>
#include <iostream>
#include <string.h>
//#include <complex.h>
#include <math.h>
#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include "coherent_beamformer.h"

using namespace std;

// Generate simulated data
float* simulate_data();

// Generate simulated weights or coefficients
float* simulate_coefficients();

// Perform transpose on the data and convert to floats
__global__
void data_transpose(float* data_in, cuComplex* data_tra);

// Convert weights from float to cuComplex
//__global__
//void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex);

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data);

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power);

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
//#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
//#endif
  return result;
}

//float* h_data = NULL;
//float* h_coeff = NULL;
float* d_data_float = NULL;
cuComplex* d_data_comp = NULL;
float* d_coeff = NULL;
//cuComplex* d_coeff_comp = NULL;
//float* d_coh_bf_out = NULL;
float* d_coh_bf_pow = NULL;
// Allocate memory to all arrays 
void init_beamformer() {
	printf("Here In init_beamformer()! \n");
	// Allocate pinned memory for input data
	//checkCuda(cudaMallocHost((void **)&h_data, N_INPUT * sizeof(float)));
	//printf("Here 1st cudaMallocHost! \n");
	
	// Allocate pinnted memery for beamformer coefficients
	//checkCuda(cudaMallocHost((void **)&h_coeff, N_COEFF * sizeof(float)));
	//printf("Here 2nd cudaMallocHost! \n");
	//cudaHostAlloc((void **)&h_data, N_INPUT * sizeof(float));

	// Allocate memory for input data float type
	checkCuda(cudaMalloc((void **)&d_data_float, N_INPUT * sizeof(float)));
	printf("Here 1st cudaMalloc! \n");

	// Allocate memory for input data cuComplex type
	checkCuda(cudaMalloc((void **)&d_data_comp, N_INPUT * sizeof(cuComplex) / 2));
	printf("Here 2nd cudaMalloc! \n");

	/*
	size_t f, t;
    	cudaSetDevice(0);
    	cudaMemGetInfo(&f, &t);
    	fprintf(stdout,"Free: %zu bytes, Available: %zu bytes \n",f,t);
	*/

	// Allocate memory for coefficients float type
	checkCuda(cudaMalloc((void **)&d_coeff, N_COEFF * sizeof(float)));
	printf("Here 3rd cudaMalloc! \n");

	// Allocate memory for coefficients cuComplex type
	//checkCuda(cudaMalloc((void **)&d_coeff_comp, N_COEFF * sizeof(cuComplex) / 2));
	//printf("Here 4th cudaMalloc! \n");

	// Allocate memory for coherent beamformer output
	//checkCuda(cudaMalloc((void **)&d_coh_bf_out, N_OUTPUT * sizeof(float)));
	//printf("Here 5th cudaMalloc! \n");

	// Allocate memory for output power of coherent beamformer
        checkCuda(cudaMalloc((void **)&d_coh_bf_pow, N_BF_POW * sizeof(float)));
	printf("Here 4th cudaMalloc! \n");

	return;
}

// Perform transpose on the data and convert to floats
__global__
void data_transpose(float* data_in, cuComplex* data_tra) {
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	data_tra[data_tr_idx(a, p, f, t)].x = data_in[2*data_in_idx(a, p, f, t)];
	data_tra[data_tr_idx(a, p, f, t)].y = data_in[2*data_in_idx(a, p, f, t) + 1];
	
	return;
}

/*
// Convert weights from float to cuComplex
__global__
void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex) {
	// Product of antenna and beam dimensions exceeds 1024 so beams
	// are blocks rather than threads to allow for increase in numbers  
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int b = blockIdx.y;  // Beam index
	int f = blockIdx.x;  // Frequency bin index 
	coeff_complex[coeff_idx(a, p, b, f)].x = coeff_float[2*coeff_idx(a, p, b, f)];
	coeff_complex[coeff_idx(a, p, b, f)].y = coeff_float[2*coeff_idx(a, p, b, f) + 1];
	
	return;
}
*/

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data) {
	/*
	int p = threadIdx.x; // Polarization index
	int f = blockIdx.x;  // Frequency index
	int t = blockIdx.y;  // Time sample index
	int b = blockIdx.z;  // Beam index
	cuComplex bf_product;
	cuComplex bf_in_data;
	cuComplex bf_coeff;
	for (int a = 0; a < N_ANT; a++) { // Antenna index
		bf_in_data.x = input_data[data_tr_idx(a, p, f, t)].x;
		bf_in_data.y = input_data[data_tr_idx(a, p, f, t)].y;
		bf_coeff.x = coeff[coeff_idx(a, p, b, f)].x;
		bf_coeff.y = coeff[coeff_idx(a, p, b, f)].y;
		// Complex multiplication of data and coefficients
		bf_product.x = (bf_in_data.x * bf_coeff.x) - (bf_in_data.y * bf_coeff.y);
		bf_product.y = (bf_in_data.x * bf_coeff.y) + (bf_in_data.y * bf_coeff.x);
		// Beamform (Sum all antennas)
		output_data[2*coh_bf_idx(p, b, f, t)] += bf_product.x;
		output_data[2*coh_bf_idx(p, b, f, t) + 1] += bf_product.y;
	}
	*/
	int a = threadIdx.x; // Antenna index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index
	int b = blockIdx.z;  // Beam index

	__shared__ cuFloatComplex reduced_mul[N_ANT];

	for (int p = 0; p < N_POL; p++) { // Polarization index
		// Reinitialize output_data since we are using the input data array to be more efficient
		int h = coh_bf_idx(p, b, f, t);
		output_data[2 * h] = 0;
		output_data[2 * h + 1] = 0;


		int i = data_tr_idx(a, p, f, t);
		int w = coeff_idx(a, b);

		if (a < N_ANT) {
			reduced_mul[a].x = input_data[i].x * coeff[2*w] + input_data[i].y * coeff[2*w + 1];
			reduced_mul[a].y = input_data[i].y * coeff[2*w] - input_data[i].x * coeff[2*w + 1];
		}
		else {
			reduced_mul[a].x = 0;
			reduced_mul[a].y = 0;
		}
		__syncthreads();

		for (int k = blockDim.x / 2; k > 0; k >>= 1) {
			if (a < k) {
				reduced_mul[a].x += reduced_mul[a + k].x;
				reduced_mul[a].y += reduced_mul[a + k].y;
			}
			__syncthreads();
		}
		if (a == 0) {
			int h = coh_bf_idx(p, b, f, t);
			output_data[2 * h] += reduced_mul[0].x;
			output_data[2 * h + 1] += reduced_mul[0].y;
		}
		
	}
	return;
}

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power) {
	int b = threadIdx.x; // Beam index
	int f = blockIdx.y;  // Frequency bin index
	int t = blockIdx.x;  // Time sample index

	// Power = Absolute value squared of output -> r^2 + i^2
	int xp = coh_bf_idx(0, b, f, t); // X polarization
	int yp = coh_bf_idx(1, b, f, t); // Y polarization
	
	float x_pol_pow = (bf_volt[2 * xp] * bf_volt[2 * xp]) + (bf_volt[2 * xp + 1] * bf_volt[2 * xp + 1]); // XX*
	float y_pol_pow = (bf_volt[2 * yp] * bf_volt[2 * yp]) + (bf_volt[2 * yp + 1] * bf_volt[2 * yp + 1]); // YY*

	bf_power[pow_bf_idx(b, f, t)] = x_pol_pow + y_pol_pow; // XX* + YY*

	/*
	bf_power[pow_bf_idx(0, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*xp]) + (bf_volt[2*xp + 1]*bf_volt[2*xp + 1]); // XX*
	bf_power[pow_bf_idx(1, b, f, t)] = (bf_volt[2*yp]*bf_volt[2*yp]) + (bf_volt[2*yp + 1]*bf_volt[2*yp + 1]); // YY*
	bf_power[pow_bf_idx(2, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*yp]) + (bf_volt[2*xp + 1]*bf_volt[2*yp + 1]); // XY* real
	bf_power[pow_bf_idx(3, b, f, t)] = (bf_volt[2*xp + 1]*bf_volt[2*yp]) - (bf_volt[2*xp]*bf_volt[2*yp + 1]); // XY* imag
	*/

	return;
}

// Run beamformer
void run_beamformer(float* data_in, float* h_coefficient, float* data_out) {
	/*
	// Allocate input data in pinned memory
	// (This may take longer than it's worth to implement pinned memory)
	*h_data = *data_in;
	*/

	cudaError_t err_code;

	// Transpose kernel: Specify grid and block dimensions
	dim3 dimBlock_transpose(N_ANT, N_POL, 1);
	dim3 dimGrid_transpose(N_TIME, N_FREQ, 1);

	// Beamformer coefficient kernel (float to complex): Specify grid and block dimensions
	//dim3 dimBlock_bf_coeff(N_ANT, N_POL, 1);
	//dim3 dimGrid_bf_coeff(N_FREQ, N_BEAM, 1);

	// Coherent beamformer kernel: Specify grid and block dimensions
	//dim3 dimBlock_coh_bf(N_ANT, N_POL, 1);
	dim3 dimBlock_coh_bf(N_ANT, 1, 1);
	dim3 dimGrid_coh_bf(N_TIME, N_FREQ, N_BEAM);

	// Output power of beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_bf_pow(N_BEAM, 1, 1);
	dim3 dimGrid_bf_pow(N_TIME, N_FREQ, 1);

	float* d_data_in = d_data_float;
	cuComplex* d_data_tra = d_data_comp;
	float* d_coefficient = d_coeff;
	//cuComplex* d_coeff_c = d_coeff_comp;
	//float* d_bf_output = d_coh_bf_out;
	float* d_bf_pow = d_coh_bf_pow;
        
	//printf("data_in[0] = %f \n", data_in[0]);
	//printf("h_coefficient[0] = %f \n", h_coefficient[0]);
	// Copy input data from host to device
	checkCuda(cudaMemcpy(d_data_in, data_in, N_INPUT * sizeof(float), cudaMemcpyHostToDevice));
	printf("First cudaMemcpy(HtoD) in run_beamformer() \n");

	// Copy beamformer coefficients from host to device
	checkCuda(cudaMemcpy(d_coefficient, h_coefficient, N_COEFF * sizeof(float), cudaMemcpyHostToDevice));
	printf("Here cudaMemcpy(HtoD) coefficients! \n");

	// Perform transpose on the data and convert to floats  
	data_transpose<<<dimGrid_transpose, dimBlock_transpose >>>(d_data_in, d_data_tra);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("BF: data_transpose() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}
	printf("Here data_transpose! \n");

	// Convert weights from float to cuComplex    
	//beamformer_coefficient<<<dimGrid_bf_coeff, dimBlock_bf_coeff >>>(d_coeff_f, d_coeff_c);
	//printf("Here beamformer_coefficient! \n");

	// Test for d_data_in
	//cudaFree(d_data_in);
	//float* d_bf_output = NULL;
	//checkCuda(cudaMalloc((void **)&d_bf_output, N_OUTPUT * sizeof(float)));

	// Perform beamforming operation
	// Use d_data_in for output since it is no longer being utilized,
	// and it is the same size as the output (4 GiB).
	coherent_beamformer<<<dimGrid_coh_bf, dimBlock_coh_bf>>>(d_data_tra, d_coefficient, d_data_in);
	//coherent_beamformer<<<dimGrid_coh_bf, dimBlock_coh_bf>>>(d_data_tra, d_coefficient, d_bf_output);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("BF: coherent_beamformer() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}
	printf("Here coherent_beamformer! \n");	

	// Compute power of beamformer output (abs()^2)    
	beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow>>>(d_data_in, d_bf_pow);
	//beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow>>>(d_bf_output, d_bf_pow);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("BF: beamformer_power() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}
	printf("Here beamformer_power! \n");

	// Copy output power from device to host
	//checkCuda(cudaMemcpy(data_out, d_bf_output, N_OUTPUT*sizeof(float), cudaMemcpyDeviceToHost));
	//checkCuda(cudaMemcpy(data_out, d_data_tra, N_INPUT*sizeof(float), cudaMemcpyDeviceToHost));
	checkCuda(cudaMemcpy(data_out, d_bf_pow, N_BF_POW*sizeof(float), cudaMemcpyDeviceToHost));
	printf("Here cudaMemcpy(DtoH)! \n");

	/*
	// Option to copy output power or voltage to host
	if(pow_flag == 0){
		cudaMemcpy(data_out, d_bf_output, N_OUTPUT*sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAerr(kern_idx);
	}else{
		beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow>>>(d_bf_output, d_bf_pow);
		checkCUDAerr(kern_idx);
		cudaMemcpy(data_out, d_bf_pow, N_BF_POW*sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAerr(kern_idx);
	}
	*/

	return;
}

// Generate simulated data
float* simulate_data() {
	float* data_sim;
	data_sim = (float*)calloc(N_INPUT, sizeof(float));
	//checkCuda(cudaMallocHost((void **)&data_sim, N_INPUT * sizeof(float)));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Repeating sequence of 1 to 64
	sim_flag = 2 -> Sequence of 1 to 64 placed in a particular bin (bin 6 for now)
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

		double theta = 0; // SOI direction/angle of arrival
		double tau = 0; // Delay
		double rf_freqs = 0;

		for (int t = 0; t < N_TIME; t++) {
			theta = ((t/100 - (N_TIME / 200)) + 90)*PI/180; // SOI direction/angle of arrival -> Moving across array over time i.e. angle changes each time sample
			tau = d * cos(theta) / c; // Delay
			for (int f = 0; f < N_FREQ; f++) {
				rf_freqs = chan_band * f + c_freq;
				for (int a = 0; a < N_ANT; a++) {
					// X polarization
					data_sim[2 * data_in_idx(a, 0, f, t)] = (float)cos(2 * PI * rf_freqs * a * tau);
					data_sim[2 * data_in_idx(a, 0, f, t) + 1] = (float)sin(2 * PI * rf_freqs * a * tau);
					// Y polarization
					data_sim[2 * data_in_idx(a, 1, f, t)] = (float)cos(2 * PI * rf_freqs * a * tau);
					data_sim[2 * data_in_idx(a, 1, f, t) + 1] = (float)sin(2 * PI * rf_freqs * a * tau); // Make this negative if a different polarization is tested
				}
			}
		}
	}
	//for (int a = 0; a < N_ANT; a++){
	//	printf("Antenna %d = %f\n", (a+1), data_sim[2 * data_in_idx(a, 0, 0, 0)]);
	//}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients() {
	float* coeff_sim;
	coeff_sim = (float*)calloc(N_COEFF, sizeof(float));
	//checkCuda(cudaMallocHost((void **)&coeff_sim, N_COEFF * sizeof(float)));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64
	sim_flag = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2
	sim flag = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.
	*/
	int sim_flag = 3;
	if (sim_flag == 0) {
		for (int i = 0; i < (N_COEFF / 2); i++) {
			coeff_sim[2 * i] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		
		for (int a = 0; a < N_ANT; a++) {
			for (int b = 0; b < N_BEAM; b++) {
				if (tmp >= N_BEAM) {
					tmp = 0;
				}
				tmp = (tmp + 1) % (N_BEAM + 1);
				coeff_sim[2 * coeff_idx(a, b)] = tmp;
			}
		}
		
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int a = 0; a < N_ANT; a++) {
			for (int b = 0; b < N_BEAM; b++) {
				if (tmp >= N_BEAM) {
					tmp = 0;
				}
				tmp = (tmp + 1) % (N_BEAM + 1);
				coeff_sim[2 * coeff_idx(a, b)] = tmp;
			}
		}
	}
	if (sim_flag == 3) {
		float c = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c / c_freq; // Wavelength
		float d = lambda / 2; // Distance between antennas
		//float chan_band = 1.59; // Fine channel bandwidth in Hz

		//float* rf_freqs = (float*)calloc(N_FREQ, sizeof(float));
		//for (int i = 0; i < N_FREQ; i++) {
		//	rf_freqs[i] = chan_band * i + c_freq;
		//}

		//float* theta = (float*)calloc(N_TIME, sizeof(float)); // Beam angle from 58 to 122 degrees
		//float* tau_beam = (float*)calloc(N_BEAM, sizeof(float)); // Delay

		double theta = 0; // Beam angle from 58 to 122 degrees
		double tau_beam = 0; // Delay

		for (int b = 0; b < N_BEAM; b++) {
			theta = ((b - (N_BEAM / 2)) + 90)*PI/180; // Beam angle from 58 to 122 degrees - Given SOI at 90 deg or moving across array, the beam with the most power is beamm 33
			tau_beam = d * cos(theta) / c; // Delay
			for (int a = 0; a < N_ANT; a++) {
				coeff_sim[2 * coeff_idx(a, b)] = (float)cos(2 * PI * c_freq * a * tau_beam);
				coeff_sim[2 * coeff_idx(a, b) + 1] = (float)sin(2 * PI * c_freq * a * tau_beam);
			}
		}
	}
	//for (int b = 0; b < N_BEAM; b++){
	//	printf("Beam %d = %f\n", (b+1), coeff_sim[2 * coeff_idx(0, b)]);
	//}

	return coeff_sim;
}

// The input_data_pin() function uses cudaHostRegister() to allocate the input host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
void input_data_pin(float * data_in_pin) {
	checkCuda(cudaHostRegister(data_in_pin, N_INPUT*sizeof(float), cudaHostRegisterPortable));
}

// The output_data_pin() function uses cudaHostRegister() to allocate the output host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
void output_data_pin(float * data_out_pin) {
	checkCuda(cudaHostRegister(data_out_pin, N_BF_POW*sizeof(float), cudaHostRegisterPortable));
}

// The coefficient_pin() function uses cudaHostRegister() to allocate the input host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
//void coefficient_pin(float * coeff_pin) {
//	checkCuda(cudaHostRegister(coeff_pin, N_COEFF*sizeof(float), cudaHostRegisterPortable));
//}

// Unregister host arrays from pinned memory
void unregister_data(void * data_unregister){
	checkCuda(cudaHostUnregister(data_unregister));
}

// Free memory
void cohbfCleanup() {
	// Free up GPU memory at the end of a program
	//if (h_data != NULL) {
	//	cudaFreeHost(h_data);
	//}
	if (d_data_float != NULL) {
		cudaFree(d_data_float);
	}
	if (d_data_comp != NULL) {
		cudaFree(d_data_comp);
	}
	if (d_coeff != NULL) {
		cudaFree(d_coeff);
	}
	//if (d_coeff_comp != NULL) {
	//	cudaFree(d_coeff_comp);
	//}
	//if (d_coh_bf_out != NULL) {
	//	cudaFree(d_coh_bf_out);
	//}
	if (d_coh_bf_pow != NULL) {
		cudaFree(d_coh_bf_pow);
	}
}


// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");

	// Allocate memory to all arrays used by run_beamformer() 
	init_beamformer();


	printf("Here1!\n");
	// Generate simulated data
	float* sim_data = simulate_data();
	// Register the array in pinned memory to speed HtoD mem copy
	input_data_pin(sim_data);

	printf("Here2!\n");

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients();
	//printf("Here3!\n");
	// Register the array in pinned memory to speed HtoD mem copy
	//coefficient_pin(sim_coefficients);

	printf("Here4!\n");

	// Allocate memory for output array
	float* output_data;
	output_data = (float*)calloc(N_BF_POW, sizeof(float));
	//output_data = (float*)calloc(N_OUTPUT, sizeof(float));
	output_data_pin(output_data);

	printf("Here5!\n");

	// Run beamformer 
	run_beamformer(sim_data, sim_coefficients, output_data);
	//run_beamformer(h_data, h_coeff, output_data);

	printf("Here6, Beamformer output: %f \n", output_data[0]);
	
	// Write data to text file for analysis
	char output_filename[128];

	printf("Here7!\n");

	//strcpy(output_filename, "C:\Users\ruzie\OneDrive\Desktop\Work\CUDA_code\output_d.txt");
	strcpy(output_filename, "output_d_cuda.txt");

	printf("Here8!\n");

	FILE* output_file;

	printf("Here9!\n");

	output_file = fopen(output_filename, "w");

	printf("Here10!\n");

	//for (int ii = 0; ii < N_INPUT; ii++) {
	for (int ii = 0; ii < N_BF_POW; ii++) {
		fprintf(output_file, "%g\n", output_data[ii]);
	}

	printf("Here11!\n");

	fclose(output_file);

	printf("Closed output file.\n");

	//unregister_data((float *)sim_data);
	//free(sim_data);
	printf("After unregister.\n");	
	free(sim_coefficients);
	printf("After freeing coefficients.\n");
	//free(output_data);

	printf("Freed dynamically allocated arrays.\n");

	// Free up device memory
	//cudaFreeHost(h_data);
	//cudaFreeHost(h_coeff);
	cohbfCleanup();

	printf("Here11!\n");

	return 0;
}
