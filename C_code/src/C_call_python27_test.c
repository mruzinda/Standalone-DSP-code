// This script tests reading and writing blocks of data to and from both RAW and text files
// Also using this script to test python module importing.
// To compile it:
// gcc C_call_python27_test.c -o C_call_python27.exe -lm -I/usr/include/python2.7 -lpython2.7
// To run it:
// ./C_call_python27.exe

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

int main()
{
	
/*
	// run objects with low-level calls //
  	int arg1 = 1;
	float *fstr;
	printf("Here \n");
  	PyObject *pmod, *pclass, *pargs, *pinst, *pmeth, *pres;
	printf("Here1 \n");

  	// instance = module.klass(  ) //
  	Py_Initialize(  );
	printf("Here2 \n");
  	pmod   = PyImport_ImportModule("func_for_C_script");         // fetch module //
	printf("Here3 \n");
  	pclass = PyObject_GetAttrString(pmod, "Vals_to_generate");   // fetch module.class //
	printf("Here4 \n");
  	Py_DECREF(pmod);
	printf("Here5 \n");

	pargs  = Py_BuildValue("(  )");
	printf("Here6 \n");
	pinst  = PyEval_CallObject(pclass, pargs);        // call class(  ) //
	printf("Here7 \n");
	Py_DECREF(pclass);
	printf("Here8 \n");
	Py_DECREF(pargs);
	printf("Here9 \n");

	// result = instance.method(x,y) //
	pmeth  = PyObject_GetAttrString(pinst, "gen_some_vals"); // fetch bound method //
	printf("Here10 \n");
	Py_DECREF(pinst);
	printf("Here11 \n");
	pargs  = Py_BuildValue("(i)", &arg1);       // convert to Python //
	printf("Here12 \n");
	//pargs  = Py_BuildValue("(ss)", arg1, arg2);       // convert to Python //
	pres   = PyEval_CallObject(pmeth, pargs);         // call method(x,y) //
	printf("Here13 \n");
	Py_DECREF(pmeth);
	printf("Here14 \n");
	Py_DECREF(pargs);
	printf("Here15 \n");
	
	PyArg_Parse(pres, "f", &fstr);                    // convert to C //
	printf("Here16 \n");
	printf("%f\n", *fstr);
	printf("Here17 \n");
	Py_DECREF(pres);
	printf("Here18 \n");
*/

	// Initialize python the python interpreter //
	Py_Initialize();

	// Set path to import module //
	PyObject *sysmodule = PyImport_ImportModule("sys");
  	assert(sysmodule != NULL);
  	PyObject *syspath = PyObject_GetAttrString(sysmodule, "path");
  	assert(syspath != NULL);
  	PyList_Append(syspath, PyString_FromString("."));
  	Py_DECREF(syspath);
  	Py_DECREF(sysmodule);

	// Import python module //
	PyObject* myModuleString = PyString_FromString((char*)"func_for_C_script");
	assert(myModuleString != NULL);
	PyObject* myModule = PyImport_Import(myModuleString);
	assert(myModule != NULL);

	// Get referrence to class //
	PyObject* myClass = PyObject_GetAttrString(myModule, (char*)"Vals_to_generate");
	assert(myClass != NULL);
	Py_DECREF(myModule);

	int arg = 1;
	int arg_flag = 0.0; // This is the argument to change the flag in the python script //
	// First argument is the size of the tuple (number of arguments).
	// Second and onward arguments are the arguments to the __init__ function of the class.
	PyObject* arglist = PyTuple_Pack(arg, PyFloat_FromDouble(arg_flag)); 
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
	Py_DECREF(myResult);

	return 0;
}
