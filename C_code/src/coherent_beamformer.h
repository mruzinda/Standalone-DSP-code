#include <stdio.h>
#include <stdlib.h>
//#include <complex.h>
#include <math.h>


#define N_POL 2 //2    // Number of polarizations
#define N_TIME 8 // 8   // Number of time samples
#define N_FREQ 10 //2^14  // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 100
#define N_ANT 64 // 64   // Number of antennas
#define N_BEAM 64 // 64  // Number of beams

// "2" for inphase and quadrature
#define N_INPUT  2*N_POL*N_TIME*N_FREQ*N_ANT       // Size of input
#define N_COEFF  2*N_ANT*N_POL*N_BEAM*N_FREQ       // Size of beamformer coefficients
#define N_OUTPUT 2*N_POL*N_BEAM*N_FREQ*N_TIME      // Size of beamformer output
#define N_BF_POW N_BEAM*N_FREQ*N_TIME    // Size of beamformer output after abs()^2

#define N_COMMANDS 2 // Number of commands run by gnuplot

#define PI 3.14159265

// p - polarization index
// t - time index
// f - frequency index
// a - antenna index
// b - beam index
#define data_in_idx(a, p, f, t)     (p + N_POL*t + N_TIME*N_POL*f + N_FREQ*N_TIME*N_POL*a)
#define data_tr_idx(a, p, f, t)     (a + N_ANT*p + N_POL*N_ANT*f + N_FREQ*N_POL*N_ANT*t)
#define coeff_idx(a, p, b, f)       (a + N_ANT*p + N_POL*N_ANT*b + N_BEAM*N_POL*N_ANT*f)
#define coh_bf_idx(p, b, f, t)      (p + N_POL*b + N_BEAM*N_POL*f + N_FREQ*N_BEAM*N_POL*t)
#define pow_bf_idx(b, f, t)         (b + N_BEAM*f + N_FREQ*N_BEAM*t)

void init_beamformer(); // Allocate memory to all arrays 
void run_beamformer(float* data_in, float* coefficient, float* data_out); // Run beamformer