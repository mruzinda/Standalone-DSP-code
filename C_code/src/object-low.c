#include <Python.h>
#include <stdio.h>

int main(  ) {

  /* run objects with low-level calls */
  char *arg1="sir", *arg2="robin", *cstr;
  PyObject *pmod, *pclass, *pargs, *pinst, *pmeth, *pres;

  /* instance = module.klass(  ) */
  Py_Initialize(  );

  PyObject *sysmodule = PyImport_ImportModule("sys");
  assert(sysmodule != NULL);
  PyObject *syspath = PyObject_GetAttrString(sysmodule, "path");
  assert(syspath != NULL);
  PyList_Append(syspath, PyString_FromString("."));
  Py_DECREF(syspath);
  Py_DECREF(sysmodule);

  pmod   = PyImport_ImportModule("module");         /* fetch module */
  assert(pmod != NULL);
  printf("Here \n");
  pclass = PyObject_GetAttrString(pmod, "klass");   /* fetch module.class */
  assert(pclass != NULL);
  printf("Here2 \n");
  Py_DECREF(pmod);
  printf("Here3 \n");

  pargs  = Py_BuildValue("");
  assert(pargs != NULL);
  printf("Here4 \n");
  pinst  = PyEval_CallObject(pclass, pargs);        /* call class(  ) */
  assert(pinst != NULL);
  printf("Here5 \n");
  Py_DECREF(pclass);
  printf("Here6 \n");
  Py_DECREF(pargs);
  printf("Here7 \n");

  /* result = instance.method(x,y) */

  pmeth  = PyObject_GetAttrString(pinst, "method"); /* fetch bound method */
  assert(pmeth != NULL);
  printf("Here8 \n");
  Py_DECREF(pinst);
  printf("Here9 \n");
  pargs  = Py_BuildValue("(ss)", arg1, arg2);       /* convert to Python */
  assert(pargs != NULL);
  printf("Here10 \n");
  pres   = PyEval_CallObject(pmeth, pargs);         /* call method(x,y) */
  assert(pres != NULL);
  printf("Here11 \n");
  Py_DECREF(pmeth);
  printf("Here12 \n");
  Py_DECREF(pargs);
  printf("Here13 \n");

  PyArg_Parse(pres, "s", &cstr);                    /* convert to C */
  printf("Here14 \n");
  printf("%s\n", cstr);
  printf("Here15 \n");
  Py_DECREF(pres);
  printf("Here16 \n");

  return 0;
}
