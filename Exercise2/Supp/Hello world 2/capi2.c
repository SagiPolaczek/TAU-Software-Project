#define PY_SSIZE_T_CLEAN
#include <Python.h>


/*
 * No arguments.
 * Note the function signature and the use of METH_NOARGS
 * Python equivalent: f1()
 */
static PyObject *f1(PyObject *self) {

    puts("f1 needs no one!");
    /* We have to return an PyObject *, so let's return the None object */
    Py_RETURN_NONE;
}


/*
 * One arguments.
 * Note the function (module and self are the common convention) signature and the use of METH_O
 * Python equivalent: f2(arg)
 * It is better to use the way we use f3() in order to not deal with reference counting
 */
static PyObject *f2(PyObject *module, PyObject *arg) {
    assert(arg);
    /* Treat arg as a borrowed reference. More on this on Moodle */
    Py_INCREF(arg);

    if (PyLong_Check(arg)) {
        long  num = PyLong_AsLong(arg); /* Missing error checking... See https://docs.python.org/3/c-api/long.html */
        printf("I like the integer number = %ld \n", num);
    }
    else puts("Can only work with integer numbers");

    /* Treat arg as a borrowed reference. More on this later. */
    Py_DECREF(arg);
    Py_RETURN_NONE;
}


/*
 * Variable Number of Arguments. You can also use it for one variable as well.
 * Note the function signature and the use of METH_VARARGS
 * Python equivalent: ֿlst = [v1,v2]; f3(lst) Or: tp = (v1,v2); f3(tp)...
 */
static PyObject* f3(PyObject *self, PyObject *args)
{
    double v1;
    int v2;
    if(!PyArg_ParseTuple(args, "(di):f3", &v1, &v2)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*
     * This builds the answer back into a python object.
     * The format unit used is almost the same as in PyArg_Parse...()
     * Here is an example where they differ - we return a dict
     */
    return Py_BuildValue("{d:i, i:d}", v1, v2, v2, v1);
}


/*
 * Variable Number of Arguments and Keyword Arguments.
 * Note the function signature and the use of METH_VARARGS | METH_KEYWORDS
 * Python equivalent: f4(*args, **keywds) or in this case f4(v1, v2='keyword v2', v3='keyword v3', v4='keyword v4')
 */
static PyObject * f4(PyObject *self, PyObject *args, PyObject *keywds) {
    int a;
    /* Initial values for optional vars. */
    const char *b = "keyword v2";
    const char *c = "keyword v3";
    const char *d = "keyword v4";
    /* The mapping of the variables, ending with a NULL sentinel */
    static char *kwlist[] = {"v1", "v2", "v3", "v4", NULL};
    /* If we don't read the optional args they keep their initial value */
    /* args, keywds are from the f4 signature. keywds will be NULL if NO keyword arguments are provided */
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "i|sss;What are you doing?", kwlist, &a, &b, &c, &d))
        return NULL;

    printf("v1 = %i\nv2 = %s\nv3 = %s\nv4 = %s\n", a, b, c, d);

    Py_RETURN_NONE;
}


/*
 * Raising exceptions and Issuing warnings.
 */
static PyObject *raise_hell(PyObject *module, PyObject *arg) {
    int error_type;
    if(!PyArg_ParseTuple(arg, "i:raise_hell", &error_type)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    switch(error_type)
    {
        case 1:
            PyErr_SetString(PyExc_NameError, "I do not like my name.");  /* Set the message text */
            return NULL; /* Raise the error */

        case 2:
            /* Message text with a format that also returns NULL */
            return PyErr_Format(PyExc_ValueError, "I can\'t stand this number: %d", error_type);

        case 3:
            /* This will also raises the error and leaves the program flow. It will NOT continue to case 4 */
            PyErr_Format(PyExc_ValueError, "4 is my lucky number, but you typed: %d", error_type);

        case 4:
            if (PyErr_Occurred()) return NULL; /* Raise the error that already occurred */
            /*
             * PyErr_WarnEx(PyObject *category, const char *message, Py_ssize_t stack_level)
             *  stack_level is a positive number giving a number of stack frames;
             *  the warning will be issued from the currently executing line of code in that stack frame.
             *  A stack_level of 1 is the function calling PyErr_WarnEx(), 2 is the function above that, and so forth.
             */

            /* By default, most types of warnings are ***only printed the first time they occur in a given location ***
               where location is defined as the combination of module and line number.  */
            else PyErr_WarnEx(PyExc_SyntaxWarning, "I warn you", 1);
            puts("I don\'t have a lucky number");
            Py_RETURN_NONE;
        default:
            return NULL; /* Raise the error */
    }
}


/*
 * Another way to creates a variable with docstring_f3 that can be used in docstrings.
 * If Python is built without docstrings, the value will be empty.
 */
PyDoc_STRVAR(docstring_f3, "Variable Number of Arguments");

static PyMethodDef _methods[] = {
    {"f1", (PyCFunction)f1, METH_NOARGS, PyDoc_STR("No input parameters")},
    {"f2", (PyCFunction)f2, METH_O, PyDoc_STR("One positional argument")},
    {"f3", (PyCFunction)f3, METH_VARARGS, docstring_f3},
    {"f4", (PyCFunction)f4, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("args and kwargs")},
    {"raise_hell", (PyCFunction)raise_hell, METH_VARARGS, PyDoc_STR("Errors raising")},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_demo2",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_capi_demo2(void)
{
    return PyModule_Create(&_moduledef);
}


