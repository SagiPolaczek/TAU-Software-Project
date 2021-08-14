#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


static PyObject* fit_spk(PyObject* self, PyObject* args);
static PyObject* fit_wam(PyObject* self, PyObject* args);
static PyObject* fit_ddg(PyObject* self, PyObject* args);
static PyObject* fit_lnorm(PyObject* self, PyObject* args);
static PyObject* fit_jacobi(PyObject* self, PyObject* args);
DataWrapper1D py_list_to_1d_array(PyObject* py_list); // TODO!
DataWrapper2D py_list_to_2d_array(PyObject* py_list);

/* Goal: SPK */
static PyObject* fit_spk(PyObject* self, PyObject* args) {
    
}

/* Goal: WAM */
static PyObject* fit_wam(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_vertices, *py_weights, *py_lnorm, *py_degrees;
    PyObject *py_result, *py_vector; // Assume WAM returns the matrix W
    double **vertices, **weights, **lnorm, *degrees, **result;
    int N, dim, K, max_iter;
    int i, j;
    DataWrapper2D vertices_wrapper, weights_wrapper, lnorm_wrapper; 
    DataWrapper1D degrees_wrapper;
    Graph G;

    if (!PyArg_ParseTuple(args, "OOOOiiii", &py_vertices, &py_weights, &py_lnorm,
    	&py_degrees, &N, &dim, &K, &max_iter))
        return NULL;
    
    /* Validation */
    assert(PyList_Check(py_vertices));
    assert(PyList_Check(py_weights));
    assert(PyList_Check(py_lnorm));
    assert(PyList_Check(py_degrees));

    /* Convert Arrays from python to C */
    vertices_wrapper = py_list_to_2d_array(py_vertices);
    weights_wrapper = py_list_to_2d_array(py_weights);
    lnorm_wrapper = py_list_to_2d_array(py_lnorm); 
    degrees_wrapper = py_list_to_1d_array(py_degrees);

    vertices = vertices_wrapper.pointers;
    weights = weights_wrapper.pointers;
    lnorm = lnorm_wrapper.pointers;
    degrees = degrees_wrapper.pointers;

    /* Initialize G */
    G.vertices = vertices;
    G.weights = weights;
    G.lnorm = lnorm;
    G.degrees = degrees;
    G.N = N;
    G.dim = dim;

    /* Main Algorithm - Assuming returns matrix W */
    /* result = compute_wam(&G); */ 
    /* Does not compile since compute_wam is void */

    /* NEED TO FINISH */

}

/* Goal: DDG */
static PyObject* fit_ddg(PyObject* self, PyObject* args) {

}

/* Goal: LNORM */
static PyObject* fit_lnorm(PyObject* self, PyObject* args) {
    
}

/* Goal: JACOBI */
static PyObject* fit_jacobi(PyObject* self, PyObject* args) {
    /* mem alloc */

    
}

/* Convert PyObject array into a double 2d array */
DataWrapper2D py_list_to_2d_array(PyObject* py_list) {
    int n, m, i, j;
    double **result;
    double *p;
    PyObject *vector, *value;
    DataWrapper2D data_wrapper;

    n = PyList_Size(py_list);
    m = PyList_Size(PyList_GetItem(py_list, 0));

    /* Init a 2-dimentaional array for the result */ 
    p = calloc((n * m), sizeof(double));
    assert(p != NULL);
    result = calloc(n, sizeof(double*));
    assert(result != NULL);

    for (i = 0; i < n; i++) {
        result[i] = p + i*m;
    }
    
    /* Put the data from the PyObject into the Array */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(py_list, i);
        assert(PyList_Check(vector));
        for (j = 0; j < m; j++) {
            value = PyList_GetItem(vector, j);
            assert(PyFloat_Check(value));
            result[i][j] = PyFloat_AsDouble(value);
        }
    }
    
    /* Pass the two array so we can free them both later */
    data_wrapper.container = p;
    data_wrapper.pointers = result;
    return data_wrapper;
}
