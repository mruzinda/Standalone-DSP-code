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
#include <device_launch_parameters.h>
#include "coherent_beamformer.h"

using namespace std;

// Check for CUDA error
//inline void checkCUDAerr(int kernel_idx);
void checkCUDAerr(int kernel_idx);

// Generate simulated data
float* simulate_data();

// Generate simulated weights or coefficients
float* simulate_coefficients();

// Perform transpose on the data and convert to floats
__global__
void data_transpose(float* data_in, cuComplex* data_tra);

// Convert weights from float to cuComplex
__global__
void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex);

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, cuComplex* coeff, float* output_data);

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power);

// Check for CUDA error
//inline void checkCUDAerr(int kernel_idx){
void checkCUDAerr(int kernel_idx) {
	//kernel_idx = kernel_idx + 1;
	cudaError_t err_code;
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		if (kernel_idx == 0) {
			printf("COH_BF: Input data cudaMemcpy() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 1) {
			printf("COH_BF: Coefficient cudaMemcpy() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 2) {
			printf("COH_BF: data_transpose() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 3) {
			printf("COH_BF: beamformer_coefficient() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 4) {
			printf("COH_BF: coherent_beamformer() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 5) {
			printf("COH_BF: beamformer_power() Failed: %s\n", cudaGetErrorString(err_code));
		}
		if (kernel_idx == 6) {
			printf("COH_BF: Final cudaMemcpy Failed: %s\n", cudaGetErrorString(err_code));
		}
	}
	return;
}

//float* h_data = NULL;
float* d_data_float = NULL;
cuComplex* d_data_comp = NULL;
float* d_coeff_float = NULL;
cuComplex* d_coeff_comp = NULL;
float* d_coh_bf_out = NULL;
float* d_coh_bf_pow = NULL;
// Allocate memory to all arrays 
void init_beamformer() {
	//cudaMallocHost((void **)&h_data, N_INPUT*sizeof(float));
	// Allocate memory for input data float type
	cudaMalloc((void**)&d_data_float, N_INPUT * sizeof(float));
	// Allocate memory for input data cuComplex type
	cudaMalloc((void**)&d_data_comp, N_INPUT * sizeof(cuComplex) / 2);
	// Allocate memory for coefficients float type
	cudaMalloc((void**)&d_coeff_float, N_COEFF * sizeof(float));
	// Allocate memory for coefficients cuComplex type
	cudaMalloc((void**)&d_coeff_comp, N_COEFF * sizeof(cuComplex) / 2);
	// Allocate memory for coherent beamformer output
	cudaMalloc((void**)&d_coh_bf_out, N_OUTPUT * sizeof(float));
	// Allocate memory for output power of coherent beamformer
	cudaMalloc((void**)&d_coh_bf_pow, N_BF_POW * sizeof(float));

	return;
}

// Perform transpose on the data and convert to floats
__global__
void data_transpose(float* data_in, cuComplex* data_tra) {
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.x;  // Frequency index
	int t = blockIdx.y;  // Time sample index

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	data_tra[data_tr_idx(a, p, f, t)].x = data_in[2 * data_in_idx(a, p, f, t)];
	data_tra[data_tr_idx(a, p, f, t)].y = data_in[2 * data_in_idx(a, p, f, t) + 1];
	
	return;
}

// Convert weights from float to cuComplex
__global__
void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex) {
	// Product of antenna and beam dimensions exceeds 1024 so beams
	// are blocks rather than threads to allow for increase in numbers  
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int b = blockIdx.y;  // Beam index
	int f = blockIdx.x;  // Frequency bin index 

	coeff_complex[coeff_idx(a, p, b, f)].x = coeff_float[2 * coeff_idx(a, p, b, f)];
	coeff_complex[coeff_idx(a, p, b, f)].y = coeff_float[2 * coeff_idx(a, p, b, f) + 1];
	
	return;
}

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, cuComplex* coeff, float* output_data) {
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
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.x;  // Frequency index
	int t = blockIdx.y;  // Time sample index
	int b = blockIdx.z;  // Beam index

	int i = data_tr_idx(a, p, f, t);
	int w = coeff_idx(a, p, b, f);

	__shared__ cuFloatComplex reduced_mul[N_ANT];

	if (a < N_ANT) {
		reduced_mul[a].x = input_data[i].x*coeff[w].x + input_data[i].y*coeff[w].y;
		reduced_mul[a].y = input_data[i].y*coeff[w].x - input_data[i].x*coeff[w].y;
	}
	else {
		reduced_mul[a].x = 0;
		reduced_mul[a].y = 0;
	}
	__syncthreads();

	for (int k = blockDim.x/2; k > 0; k >>= 1) {
		if (a < k) {
			reduced_mul[a].x += reduced_mul[a + k].x;
			reduced_mul[a].y += reduced_mul[a + k].y;
		}
		__syncthreads();
	}
	if (a == 0) {
		output_data[2*coh_bf_idx(p, b, f, t)] += reduced_mul[0].x;
		output_data[2*coh_bf_idx(p, b, f, t) + 1] += reduced_mul[0].y;
	}

	return;
}

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power) {
	int b = threadIdx.x; // Beam index
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.x;  // Frequency bin index
	int t = blockIdx.y;  // Time sample index

	int h = coh_bf_idx(p, b, f, t);
	// Power = Absolute value squared of output -> r^2 + i^2
	bf_power[h] = (bf_volt[2*h]*bf_volt[2*h]) + (bf_volt[2*h + 1]*bf_volt[2*h + 1]);
	
	return;
}

// Run beamformer
void run_beamformer(float* data_in, float* h_coefficient, float* data_out) {
	int kern_idx = 0; // Kernel index in run_beamformer function for printing error

	/*
	// Allocate input data in pinned memory
	// (This may take longer than it's worth to implement pinned memory)
	*h_data = *data_in;
	*/

	// Transpose kernel: Specify grid and block dimensions
	dim3 dimBlock_transpose(N_ANT, N_POL, 1);
	dim3 dimGrid_transpose(N_BIN, N_TIME, 1);

	// Beamformer coefficient kernel (float to complex): Specify grid and block dimensions
	dim3 dimBlock_bf_coeff(N_ANT, N_POL, 1);
	dim3 dimGrid_bf_coeff(N_BIN, N_BEAM, 1);

	// Coherent beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_coh_bf(N_ANT, N_POL, 1);
	dim3 dimGrid_coh_bf(N_BIN, N_TIME, N_BEAM);

	// Output power of beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_bf_pow(N_BEAM, N_POL, 1);
	dim3 dimGrid_bf_pow(N_BIN, N_TIME, 1);

	float* d_data_in = d_data_float;
	cuComplex* d_data_tra = d_data_comp;
	float* d_coeff_f = d_coeff_float;
	cuComplex* d_coeff_c = d_coeff_comp;
	float* d_bf_output = d_coh_bf_out;
	float* d_bf_pow = d_coh_bf_pow;

	// Copy input data from host to device
	cudaMemcpy(d_data_in, data_in, N_INPUT * sizeof(float), cudaMemcpyHostToDevice);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Copy beamformer coefficients from host to device
	cudaMemcpy(d_coeff_f, h_coefficient, N_COEFF * sizeof(float), cudaMemcpyHostToDevice);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Perform transpose on the data and convert to floats  
	data_transpose<<<dimGrid_transpose, dimBlock_transpose >>>(d_data_in, d_data_tra);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Convert weights from float to cuComplex    
	beamformer_coefficient<<<dimGrid_bf_coeff, dimBlock_bf_coeff >>>(d_coeff_f, d_coeff_c);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Perform beamforming operation
	coherent_beamformer<<<dimGrid_coh_bf, dimBlock_coh_bf>>>(d_data_tra, d_coeff_c, d_bf_output);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Compute power of beamformer output (abs()^2)    
	beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow>>>(d_bf_output, d_bf_pow);
	checkCUDAerr(kern_idx);
	kern_idx = kern_idx + 1;
	printf("Here kernel idx: %d!\n", kern_idx);

	// Copy output power from device to host
	cudaMemcpy(data_out, d_bf_pow, N_BF_POW*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAerr(kern_idx);
	printf("Here kernel idx: %d!\n", kern_idx);

	kern_idx = 0; // Reset the kernel index for the CUDA error check
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
	for (int i = 0; i < (N_INPUT / 2); i++) {
		data_sim[2*i] = 1;
	}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients() {
	float* coeff_sim;

	printf("Here in sim coeff!\n");

	coeff_sim = (float*)calloc(N_COEFF, sizeof(float));

	printf("Here in sim coeff 1!\n");

	for (int i = 0; i < (N_COEFF / 2); i++) {
		coeff_sim[2*i] = 1;
		//if(i> (3*N_COEFF/4)){
		//printf("Here in sim coeff: %f, idx = %d\n", coeff_sim[2*i], i);
		//}
	}

	printf("Here in sim coeff 2!\n");

	return coeff_sim;
}

// Free memory
void cohbfCleanup() {
	// Free up GPU memory at the end of a program
	if (d_data_float != NULL) {
		cudaFree(d_data_float);
	}
	if (d_data_comp != NULL) {
		cudaFree(d_data_comp);
	}
	if (d_coeff_float != NULL) {
		cudaFree(d_coeff_float);
	}
	if (d_coeff_comp != NULL) {
		cudaFree(d_coeff_comp);
	}
	if (d_coh_bf_out != NULL) {
		cudaFree(d_coh_bf_out);
	}
	if (d_coh_bf_pow != NULL) {
		cudaFree(d_coh_bf_pow);
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

	// Write data to text file for analysis
	char output_filename[128];

	printf("Here6!\n");

	//strcpy(output_filename, "C:\Users\ruzie\OneDrive\Desktop\Work\CUDA_code\output_d.txt");
	strcpy(output_filename, "output_d.txt");

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