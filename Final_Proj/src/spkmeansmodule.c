#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "spkmeans.h"
#include "kmeans.h"

#if PY_MAJOR_VERSION >= 3
#define PY3K
#endif

static void fit_general(PyObject* self, PyObject* args); /* wam ddg lnorm jacobi */
static PyObject* fit_init_spk(PyObject* self, PyObject* args);
static void fit_finish_spk(PyObject* self, PyObject* args);
void py_list_to_array(PyObject* py_list, int n, int m, DataWrapper data_points_wrapper);

static void fit_general(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    double **data_points;
    int N, dim, K, max_iter, g;
    int n, m, i;
    goal goal;
    DataWrapper data_points_wrapper;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiiii", &py_data_points, &N, &dim, &K, &max_iter, &g))
        return;
    
    /* Validation */
    assert(PyList_Check(py_data_points));

    n = (int)PyList_Size(py_data_points);
    m = (int)PyList_Size(PyList_GetItem(py_data_points, 0));

    /* Init a 2-dimentaional array for the result */ 
    data_points_wrapper.container = calloc((n * m), sizeof(double));
    assert(data_points_wrapper.container != NULL);
    data_points_wrapper.pointers = calloc(n, sizeof(double*));
    assert(data_points_wrapper.pointers != NULL);

    for (i = 0; i < n; i++) {
        data_points_wrapper.pointers[i] = data_points_wrapper.container + i*m;
    }

    /* Convert Arrays from Python to C */
    py_list_to_array(py_data_points, n, m, data_points_wrapper);
    data_points = data_points_wrapper.pointers;

    /* Init graph */
    graph.vertices = data_points;
    graph.N = N;
    graph.dim = dim;

    goal = g;
    /* Main algorithm */
    compute_by_goal(&graph, goal);

    /* Free memory */
    free(data_points_wrapper.container);
    free(data_points_wrapper.pointers);
}

static PyObject* fit_init_spk(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    PyObject *py_result, *py_vector;
    double **data_points;
    double **result;
    int N, dim, K, max_iter;
    int i, j, n, m;
    DataWrapper data_points_wrapper;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiii", &py_data_points, &N, &dim, &K, &max_iter))
        return NULL;
    
    /* Validation */
    assert(PyList_Check(py_data_points));

    n = (int)PyList_Size(py_data_points);
    m = (int)PyList_Size(PyList_GetItem(py_data_points, 0));

    /* Init a 2-dimentaional array for the result */ 
    data_points_wrapper.container = calloc((n * m), sizeof(double));
    assert(data_points_wrapper.container != NULL);
    data_points_wrapper.pointers = calloc(n, sizeof(double*));
    assert(data_points_wrapper.pointers != NULL);

    for (i = 0; i < n; i++) {
        data_points_wrapper.pointers[i] = data_points_wrapper.container + i*m;
    }

    /* Convert Arrays from Python to C */
    py_list_to_array(py_data_points, n, m, data_points_wrapper);
    data_points = data_points_wrapper.pointers;

    /* Init graph */
    graph.vertices = data_points;
    graph.N = N;
    graph.dim = dim;

    result = init_spk_datapoints(&graph, &K);

    /* Convert a two double array into a PyObject */
    py_result = PyList_New(K);
    assert(py_result != NULL && "Problem in generating PyList object");
    for (i = 0; i < K; i++) {
        py_vector = PyList_New(dim);
        assert(py_vector != NULL && "Problem in generating PyList object");
        for (j = 0; j < dim; j++) {
            PyList_SET_ITEM(py_vector, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SET_ITEM(py_result, i, py_vector);
    }

    free(data_points_wrapper.container);
    free(result);

    return py_result;
}

static void fit_finish_spk(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_centroids, *py_data, *py_indices;
    double **centroids, **data_points, **result;
    int N, dim, K, max_iter;
    int n, m, i;
    DataWrapper centroids_wrapper, data_points_wrapper;

    if (!PyArg_ParseTuple(args, "OOOiiii", &py_centroids, &py_data, &py_indices, &N, &dim, &K, &max_iter))
        return;
    
    /* Validation */
    assert(PyList_Check(py_centroids));
    assert(PyList_Check(py_data));
    assert(PyList_Check(py_indices));

    /* Allocate memory for data_points_wrapper */

    n = (int)PyList_Size(py_data);
    m = (int)PyList_Size(PyList_GetItem(py_data, 0));

    /* Init a 2-dimentaional array for the result */ 
    data_points_wrapper.container = calloc((n * m), sizeof(double));
    assert(data_points_wrapper.container != NULL);
    data_points_wrapper.pointers = calloc(n, sizeof(double*));
    assert(data_points_wrapper.pointers != NULL);

    for (i = 0; i < n; i++) {
        data_points_wrapper.pointers[i] = data_points_wrapper.container + i*m;
    }

    /* Convert Arrays from Python to C */
    py_list_to_array(py_data, n, m, data_points_wrapper);
    data_points = data_points_wrapper.pointers;

    /* Allocate memory for centroids_wrapper */

    n = (int)PyList_Size(py_centroids);
    m = (int)PyList_Size(PyList_GetItem(py_centroids, 0));

    /* Init a 2-dimentaional array for the result */ 
    centroids_wrapper.container = calloc((n * m), sizeof(double));
    assert(centroids_wrapper.container != NULL);
    centroids_wrapper.pointers = calloc(n, sizeof(double*));
    assert(centroids_wrapper.pointers != NULL);

    for (i = 0; i < n; i++) {
        centroids_wrapper.pointers[i] = centroids_wrapper.container + i*m;
    }

    /* Convert Arrays from Python to C */
    py_list_to_array(py_centroids, n, m, centroids_wrapper);
    centroids = centroids_wrapper.pointers;

    /* Main Algorithm - should call another func here which will call kmeans */
    result = kmeans(data_points, centroids, N, dim, K, max_iter);
    /* with the new interface:
    result = get_spk_clusters(data_points, centroids, N, dim, K, max_iter); */

    /* Print indices */
    for (i = 0; i < K; i++) {
        printf("%.4f", PyFloat_AsDouble(PyList_GET_ITEM(py_indices, i)));
        if (i < K -1) {
            printf(",");
        }
        else {
            printf("/n");
        }
    }
    /* Print centroids */
    print_matrix(result, N, K);

    /* Free Memory */
    free(centroids_wrapper.container);
    free(data_points_wrapper.container);
}

/* Convert PyObject array into a double 2d array */
void py_list_to_array(PyObject* py_list, int n, int m, DataWrapper data_points_wrapper) {
    int i, j;
    PyObject *vector, *value;
    
    /* Put the data from the PyObject into the Array */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(py_list, i);
        assert(PyList_Check(vector));
        for (j = 0; j < m; j++) {
            value = PyList_GetItem(vector, j);
            assert(PyFloat_Check(value));
            data_points_wrapper.pointers[i][j] = PyFloat_AsDouble(value);
        }
    }
}

/* Python Staff */
static PyMethodDef _methods[] = {
    {"fit_general", (PyCFunction)fit_general, METH_VARARGS, PyDoc_STR("Program for wam, ddg, lnorm and jacobi goals")},
    {"fit_init_spk", (PyCFunction)fit_init_spk, METH_VARARGS, PyDoc_STR("Initials the T matrix for spk goal")},
    {"fit_finish_spk", (PyCFunction)fit_finish_spk, METH_VARARGS, PyDoc_STR("Preforms kmeans algorithm for spk goal")},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

#ifdef PY3K
static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "myspkmeans",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC PyInit_myspkmeans(void)
{
    return PyModule_Create(&_moduledef);
}

#else
PyMODINIT_FUNC initmyspkmeans(void) {
    Py_InitModule3("myspkmeans", _methods, NULL);
}

#endif
