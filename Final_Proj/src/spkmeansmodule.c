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
void py_list_to_array(PyObject* py_list, int n, int m, double **data_points);

static void fit_general(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    double **data_points;
    int N, dim, max_iter, g;
    int n, m;
    goal goal;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiii", &py_data_points, &N, &dim, &max_iter, &g))
        return;
    printf("C-%d\n", 1);
    /* Validation */
    assert(PyList_Check(py_data_points));

    n = (int)PyList_Size(py_data_points);
    m = (int)PyList_Size(PyList_GetItem(py_data_points, 0));

    data_points = calloc_2d_array(n, m);

    /* Convert Arrays from Python to C */
    py_list_to_array(py_data_points, n, m, data_points);
    
    printf("C-%d\n", 3);
    /* Init graph */
    graph.vertices = data_points;
    graph.N = N;
    graph.dim = dim;

    goal = g;
    /* Main algorithm */
    compute_by_goal(&graph, goal);
    printf("C-%d\n", 4);
    /* Free memory */
    /* free(&graph); */
}

static PyObject* fit_init_spk(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    PyObject *py_result, *py_vector;
    double **data_points;
    double **result;
    int N, dim, K, max_iter;
    int i, j, n, m;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiii", &py_data_points, &N, &dim, &K, &max_iter))
        return NULL;
    
    /* Validation */
    assert(PyList_Check(py_data_points));

    n = (int)PyList_Size(py_data_points);
    m = (int)PyList_Size(PyList_GetItem(py_data_points, 0));

    /* Init a 2-dimentaional array */ 
    data_points = calloc_2d_array(n, m);

    /* Convert Arrays from Python to C */
    py_list_to_array(py_data_points, n, m, data_points);

    /* Init graph */
    graph.vertices = data_points;
    graph.N = N;
    graph.dim = dim;

    result = init_spk_datapoints(&graph, &K);

    /* Convert a two double array into a PyObject */
    py_result = PyList_New(N);
    assert(py_result != NULL && "Problem in generating PyList object");
    for (i = 0; i < N; i++) {
        py_vector = PyList_New(K);
        assert(py_vector != NULL && "Problem in generating PyList object");
        for (j = 0; j < K; j++) {
            PyList_SET_ITEM(py_vector, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SET_ITEM(py_result, i, py_vector);
    }

    free_2d_array(data_points);
    free(result);

    return py_result;
}

static void fit_finish_spk(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_centroids, *py_data, *py_indices;
    double **centroids, **data_points, **result;
    int N, dim, K, max_iter;
    int n, m;

    if (!PyArg_ParseTuple(args, "OOOiiii", &py_centroids, &py_data, &py_indices, &N, &dim, &K, &max_iter))
        return;
    
    /* Validation */
    assert(PyList_Check(py_centroids));
    assert(PyList_Check(py_data));
    assert(PyList_Check(py_indices));

    /* Allocate memory for data_points_wrapper */

    n = (int)PyList_Size(py_data);
    m = (int)PyList_Size(PyList_GetItem(py_data, 0));

    /* Init a 2-dimentaional array for data_points */ 
    data_points = calloc_2d_array(n, m);
    
    /* Convert Arrays from Python to C */
    py_list_to_array(py_data, n, m, data_points);

    /* Allocate memory for centroids_wrapper */

    n = (int)PyList_Size(py_centroids);
    m = (int)PyList_Size(PyList_GetItem(py_centroids, 0));

    /* Init a 2-dimentaional array for data_points */ 
    centroids = calloc_2d_array(n, m);
    
    /* Convert Arrays from Python to C */
    py_list_to_array(py_data, n, m, centroids);

    result = centroids;
    /* Main Algorithm - should call another func here which will call kmeans */
    kmeans(data_points, N, dim, K, max_iter, result);

    /* Print centroids */
    print_matrix(result, N, K);

    /* Free Memory */
    free_2d_array(data_points);
    free_2d_array(centroids);
}

/* Convert PyObject array into a double 2d array */
void py_list_to_array(PyObject* py_list, int n, int m, double **data_points) {
    int i, j;
    PyObject *vector, *value;
    
    /* Put the data from the PyObject into the Array */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(py_list, i);
        assert(PyList_Check(vector));
        for (j = 0; j < m; j++) {
            value = PyList_GetItem(vector, j);
            assert(PyFloat_Check(value));
            data_points[i][j] = PyFloat_AsDouble(value);
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
