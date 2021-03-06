// Using this script to test python module importing. Specifically, katpoint and mosaic modules for delay calculation.
// Use python 3.7
// To compile it:
// gcc Delay_katpoint_mosaic.c -o Delay_katpoint_mosaic.exe -lm -I/opt/conda/include/python3.7m -L/opt/conda/lib -lpython3.7m
// And make sure that /opt/conda/lib is include in the LD_LIBRARY_PATH
// Also have some version of these paths to your PYTHONPATH i.e. just make sure you have the paths to the modules in the PYTHONPATH:
// /home/mruzinda/Calculate_delay/Beamforming/mosaic:/home/mruzinda/Calculate_delay:/home/mruzinda/Standalone-DSP-code/Python_code:/home/mruzinda/Standalone-DSP-code/C_code/src
// To run it:
// ./Delay_katpoint_mosaic.exe

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdint.h>
#include <endian.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#define N_ANT 64
#define DELAY_POLYS 2
#define delay_idx(d, a, b)  (d + DELAY_POLYS*a + DELAY_POLYS*N_ANT*b) // Should be correct indexing

int main()
{
	float time_taken = 0;
	float delaycalc_time = 0;
	int arr_size = 0;
	float* result;
	int num_runs = 10;

	PyObject* myModuleString;
	PyObject* myModule;
	PyObject* myClass;
	PyObject* arglist;
	PyObject* myInst;
	PyObject* myMethod;
	PyObject* arglist2;
	PyObject* myResult;
	PyObject* result_tmp;

	// Start timing delay calculation //
	struct timespec tval_before, tval_after;

	for(int ii = 0; ii < num_runs; ii++){

		// Start timing delay calculation //
		clock_gettime(CLOCK_MONOTONIC, &tval_before);

		// Initialize python the python interpreter //
		Py_Initialize();
	
		// Import python module //
		myModuleString = PyUnicode_DecodeFSDefault("get_delays");
		//PyObject* myModuleString = PyUnicode_DecodeFSDefaultAndSize("test_module",11);
		//PyObject* myModuleString = PyUnicode_DecodeFSDefault("func_for_C_script");
		//PyObject* myModuleString = PyUnicode_FromString("func_for_C_script");
		assert(myModuleString != NULL);
		myModule = PyImport_Import(myModuleString);
		//PyObject* myModule = PyImport_ImportModule("func_for_C_script");
		//PyObject* myModule = PyImport_ImportModule("test_module");
		assert(myModule != NULL);
		Py_DECREF(myModuleString);
	
		// Get referrence to class //
		myClass = PyObject_GetAttrString(myModule, "DelayPolynomial");
		assert(myClass != NULL);
		Py_DECREF(myModule);

		int arg = 1;
		float freq_chan = 1.4e9; // This is the argument to change the flag in the python script //
		// First argument is the size of the tuple (number of arguments).
		// Second and onward arguments are the arguments to the __init__ function of the class.
		arglist = PyTuple_Pack(arg, PyFloat_FromDouble(freq_chan)); 
		//PyObject* arglist = PyTuple_Pack(arg, PyUnicode_DecodeFSDefault("1")); 
		assert(arglist != NULL);

		// Get class //
		myInst = PyObject_CallObject(myClass, arglist);
		assert(myInst != NULL);
	
		// Get referrence to method/function //
		myMethod  = PyObject_GetAttrString(myInst, "get_delay_polynomials"); // fetch bound method //
  		assert(myMethod != NULL);
  		Py_DECREF(myInst);

		int arg2 = 2;
		int time_arg = 1629380016;
		int dur_flag = 2;
		// No argument needed for the gen_some_vals //
		arglist2 = PyTuple_Pack(arg2, PyFloat_FromDouble(time_arg), PyFloat_FromDouble(dur_flag));
		assert(arglist2 != NULL);

		// Get result from function //
		myResult = PyObject_CallObject(myMethod, arglist2);
		Py_DECREF(myMethod);
		Py_DECREF(arglist2);
		assert(myResult != NULL);
	
		/*
		// Gets a single value from a variable in python //
		float result = (float)PyFloat_AsDouble(myResult);
		Py_DECREF(myResult);
		printf("Result from python module = %f \n", result);
		*/

		// Gets the size of the array generated by python script //
		int arr_size = PyList_Size(myResult);
		//printf("Length of array = %d \n", arr_size);

		// Gets the values from the array //
		result = (float*)calloc(arr_size, sizeof(float));
		for(int i = 0; i < arr_size; i++){
			result_tmp = PyList_GetItem(myResult, i);
			result[i] = (float)PyFloat_AsDouble(result_tmp);
		}

		// Stop timing beamforming computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
		delaycalc_time += time_taken;
		printf("Time taken: %f s\n", time_taken);
	}

	printf("Average delay calculation time: %f s\n", delaycalc_time/num_runs);

	// First beam
	printf("--------------First beam delay offset---------------\n");
	printf("idx %d in result array = %e \n", delay_idx(0, 0, 0), result[delay_idx(0, 0, 0)]);
	printf("idx %d in result array = %e \n", delay_idx(0, 1, 0), result[delay_idx(1, 0, 0)]);
	printf("idx %d in result array = %e \n", delay_idx(0, 2, 0), result[delay_idx(2, 0, 0)]);
	// Second beam delay
	printf("--------------Second beam delay offset--------------\n");
	printf("idx %d in result array = %e \n", delay_idx(0, 0, 1), result[delay_idx(0, 0, 1)]);
	printf("idx %d in result array = %e \n", delay_idx(0, 1, 1), result[delay_idx(0, 1, 1)]);
	printf("idx %d in result array = %e \n", delay_idx(0, 2, 1), result[delay_idx(0, 2, 1)]);
	// Second beam rate
	printf("---------------Second beam delay rate----------------\n");
	printf("idx %d in result array = %e \n", delay_idx(1, 0, 1), result[delay_idx(1, 0, 1)]); // 129
	printf("idx %d in result array = %e \n", delay_idx(1, 1, 1), result[delay_idx(1, 1, 1)]); // 131
	printf("idx %d in result array = %e \n", delay_idx(1, 2, 1), result[delay_idx(1, 2, 1)]); // 133

	Py_Finalize();

	return 0;
}
