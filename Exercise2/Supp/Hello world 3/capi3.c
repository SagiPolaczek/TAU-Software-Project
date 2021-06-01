#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*
 * Helper functions prototypes
 */

static void print_list_of_ints(PyObject *, Py_ssize_t);

/*
 * Helper functions implementation
 */


/*
 * Print a list of ints 
 * All integers are implemented as “long” integer objects of arbitrary size.
 */
static void print_list_of_ints(PyObject * list, Py_ssize_t list_size){
    Py_ssize_t i;
    PyObject *item;
    /* NEVER EVER USE malloc/calloc/realloc or free on PyObject */
    /* Here we use it correctly on our own C type */
    /* All integers are implemented as “long” integer objects of arbitrary size */
    long *my_c_list = malloc(sizeof my_c_list * list_size);
    /* Check that we got the memory from OS. In the assert - a string has a true value */
    assert(my_c_list != NULL && "Problem in print_list_of_ints()");
    puts("*****************************************");
    for (i = 0; i < list_size; i++) {
        item = PyList_GetItem(list, i); /* Return value: Borrowed reference - more on this in Moodle */
        if (!PyLong_Check(item)) continue; /* Skip non-integers */
        /* See here https://docs.python.org/3.8/c-api/concrete.html under Numeric Objects */
        my_c_list[i] = PyLong_AsLong(item); /* Convert a Python int object to long */
        if (my_c_list[i]  == -1 && PyErr_Occurred()){
            /* Integer too big to fit in a C long, bail out */
            puts("Something bad ...");
            /* Don't forget to free back the memory */
            free(my_c_list);
            return;
         }
        printf("%ld\n", my_c_list[i]);
    }
    puts("*****************************************");
    /* Don't forget to free back the memory */
    free(my_c_list);
    return;
}




/*
 * API functions
 */


/*
 * Print list of lists of ints without changing it
 */
static PyObject* print_int_lists(PyObject *self, PyObject *args)
{
    PyObject *_list, *item;
    Py_ssize_t i, n;

    if(!PyArg_ParseTuple(args, "O:print_int_lists", &_list)) {
        return NULL;
    }
    /* Is it a list? */
    if (!PyList_Check(_list))
        return NULL;
    /* Get the size of it and build the output list */
    n = PyList_Size(_list);  /*  Same as in Python len(_list)  */
    /* Go over each item of the list and reduce it */
    for (i = 0; i < n; i++) {
        item = PyList_GetItem(_list, i);
        if (!PyList_Check(item)){  /* We only print lists */
           continue;
        }
        print_list_of_ints(item, PyList_Size(item));
    }
    Py_RETURN_NONE;
}

/*
 * A macro to help us with defining the methods
 * Compare with: {"f1", (PyCFunction)f1, METH_NOARGS, PyDoc_STR("No input parameters")}
*/
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
    FUNC(METH_VARARGS, print_int_lists, "Print list of lists of ints without changing it"),
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_demo3",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_capi_demo3(void)
{
    return PyModule_Create(&_moduledef);
}


