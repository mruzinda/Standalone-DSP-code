#include <stdio.h>
#include <stdlib.h>
//#include <complex.h>
#include <math.h>


#define N_POL 2 //2                     // Number of polarizations
#define N_TIME 8 // 8                   // Number of time samples
#define N_COARSE_FREQ 32               // Number of coarse channels processed at a time
#define N_FINE_FREQ 16384               // Number of fine channels per coarse channel 2^14 = 16384
#define N_FREQ N_COARSE_FREQ*N_FINE_FREQ // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_ANT 64 // 64                  // Number of antennas
#define N_BEAM 64 // 64                 // Number of beams
//#define N_POL_OUT 4 //2    // Number of output polarizations 

// "2" for inphase and quadrature
#define N_INPUT  (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_ANT)       // Size of input. Currently, same size as output
//#define N_COEFF  (unsigned long int)(2*N_ANT*N_POL*N_BEAM*N_FREQ)     // Size of beamformer coefficients
#define N_COEFF  (unsigned long int)(2*N_ANT*N_BEAM)                    // Size of beamformer coefficients
#define N_OUTPUT (unsigned long int)(2*N_POL*N_BEAM*N_FREQ*N_TIME)      // Size of beamformer output
#define N_BF_POW (unsigned long int)(N_BEAM*N_FREQ*N_TIME)              // Size of beamformer output after abs()^2
//#define N_BF_POW N_POL_OUT*N_BEAM*N_FREQ*N_TIME    // Size of beamformer output after abs()^2

#define PI 3.14159265

// p - polarization index
// t - time index
// f - frequency index
// a - antenna index
// b - beam index
#define data_in_idx(a, p, f, t)     (p + N_POL*t + N_TIME*N_POL*f + N_FREQ*N_TIME*N_POL*a)
#define data_tr_idx(a, p, f, t)     (a + N_ANT*p + N_POL*N_ANT*f + N_FREQ*N_POL*N_ANT*t)
#define coeff_idx(a, b)             (a + N_ANT*b)
#define coh_bf_idx(p, b, f, t)      (p + N_POL*b + N_BEAM*N_POL*f + N_FREQ*N_BEAM*N_POL*t)
#define pow_bf_idx(b, f, t)         (b + N_BEAM*f + N_FREQ*N_BEAM*t)

/*
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
//#define pow_bf_idx(p, b, f, t)      (p + N_POL_OUT*b + N_BEAM*N_POL_OUT*f + N_FREQ*N_BEAM*N_POL_OUT*t)
*/

#ifdef __cplusplus
extern "C" {
#endif
void init_beamformer(); // Allocate memory to all arrays 
float* simulate_data();
float* simulate_coefficients();
void input_data_pin(float * data_pin);
//void coefficient_pin(float * coeff_pin);
void unregister_data(void * data_unregister);
void cohbfCleanup();
void run_beamformer(float* data_in, float* coefficient, float* data_out); // Run beamformer
#ifdef __cplusplus
}
#endif
