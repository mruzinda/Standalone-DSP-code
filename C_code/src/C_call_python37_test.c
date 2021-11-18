// Using this script to test python module importing.
// To compile it:
// gcc C_call_python35_test.c -o C_call_python35.exe -lm -I/usr/include/python3.5 -lpython3.5m
// To run it:
// ./C_call_python35.exe
// For python 3.7
// To compile it:
// gcc C_call_python37_test.c -o C_call_python37.exe -lm -I/opt/conda/include/python3.7m -L/opt/conda/lib -lpython3.7m
// And make sure that /opt/conda/lib is include in the LD_LIBRARY_PATH

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <endian.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
//#define SYSPATH L"/home/mruzinda/Standalone-DSP-code/C_code/src"
#define SYSPATH L"/home/mruzinda/Standalone-DSP-code/Python_code"

int main()
{

	// Initialize python the python interpreter //
	Py_Initialize();

	// Import python module //
	PyObject* myModuleString = PyUnicode_DecodeFSDefault("test_module");
	//PyObject* myModuleString = PyUnicode_DecodeFSDefaultAndSize("test_module",11);
	//PyObject* myModuleString = PyUnicode_DecodeFSDefault("func_for_C_script");
	//PyObject* myModuleString = PyUnicode_FromString("func_for_C_script");
	assert(myModuleString != NULL);
	PyObject* myModule = PyImport_Import(myModuleString);
	//PyObject* myModule = PyImport_ImportModule("func_for_C_script");
	//PyObject* myModule = PyImport_ImportModule("test_module");
	assert(myModule != NULL);
	Py_DECREF(myModuleString);
	
	// Get referrence to class //
	PyObject* myClass = PyObject_GetAttrString(myModule, "Vals_to_generate");
	assert(myClass != NULL);
	Py_DECREF(myModule);

	int arg = 1;
	//int arg_flag = 0.0; // This is the argument to change the flag in the python script //
	// First argument is the size of the tuple (number of arguments).
	// Second and onward arguments are the arguments to the __init__ function of the class.
	//PyObject* arglist = PyTuple_Pack(arg, PyFloat_FromDouble(arg_flag)); 
	PyObject* arglist = PyTuple_Pack(arg, PyUnicode_DecodeFSDefault("1")); 
	assert(arglist != NULL);

	// Get class //
	PyObject* myInst = PyObject_CallObject(myClass, arglist);
	assert(myInst != NULL);

	// Get referrence to method/function //
	PyObject* myMethod  = PyObject_GetAttrString(myInst, "gen_some_vals"); // fetch bound method //
  	assert(myMethod != NULL);
  	Py_DECREF(myInst);

	int arg2 = 0;
	// No argument needed for the gen_some_vals //
	PyObject* arglist2 = PyTuple_Pack(arg2, PyFloat_FromDouble(1.0));
	assert(arglist2 != NULL);

	// Get result from function //
	PyObject* myResult = PyObject_CallObject(myMethod, arglist2);
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
	printf("Length of array = %d \n", arr_size);

	// Gets the values from the array //
	PyObject* result_tmp;
	float* result = (float*)calloc(arr_size, sizeof(float));
	for(int i = 0; i < arr_size; i++){
		result_tmp = PyList_GetItem(myResult, i);
		result[i] = (float)PyFloat_AsDouble(result_tmp);
		printf("idx %d in result array = %f \n", i, result[i]);
	}

	Py_Finalize();

	return 0;
}
