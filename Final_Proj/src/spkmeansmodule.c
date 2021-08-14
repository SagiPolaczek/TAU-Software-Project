#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "spkmeans.h"
#include "kmeans.h"

static void fit_general(PyObject* self, PyObject* args); /* wam ddg lnorm jacobi */
static PyObject* fit_init_spk(PyObject* self, PyObject* args);
static void fit_finish_spk(PyObject* self, PyObject* args);
DataWrapper py_list_to_array(PyObject* py_list);

static void fit_general(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    double **data_points;
    double *eigenvalues, **eigenvectors, **A;
    int N, dim, K, max_iter, g;
    int i, j;
    goal goal;
    DataWrapper data_points_wrapper;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiiii", &py_data_points, &N, &dim, &K, &max_iter, &g))
        return;
    
    /* Validation */
    assert(PyList_Check(py_data_points));

    /* Convert Arrays from Python to C */
    data_points_wrapper = py_list_to_array(py_data_points);
    data_points = data_points_wrapper.pointers;

    /* Init graph */
    graph.vertices = data_points;
    graph.N = N;
    graph.dim = dim;

    goal = g;
    
    /* compute_result(graph, goal); */

    /* with the new interface:
    compute_by_goal(&graph, goal); */

    /* Free memory */
    free(data_points_wrapper.container);
}

static PyObject* fit_init_spk(PyObject* self, PyObject* args) {
    /* Variables Declarations */
    PyObject *py_data_points;
    PyObject *py_result, *py_vector;
    double **data_points;
    double *eigenvalues, **eigenvectors, **A, **result;
    int N, dim, K, max_iter;
    int i, j;
    DataWrapper data_points_wrapper;
    Graph graph = {0};

    if (!PyArg_ParseTuple(args, "Oiiii", &py_data_points, &N, &dim, &K, &max_iter))
        return NULL;
    
    /* Validation */
    assert(PyList_Check(py_data_points));

    /* Convert Arrays from Python to C */
    data_points_wrapper = py_list_to_array(py_data_points);
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
    PyObject *py_result, *py_vector;
    double **centroids, **data_points, **result;
    int N, dim, K, max_iter;
    int i, j;
    DataWrapper centroids_wrapper, data_points_wrapper;

    if (!PyArg_ParseTuple(args, "OOOiiii", &py_centroids, &py_data, &py_indices, &N, &dim, &K, &max_iter))
        return;
    
    /* Validation */
    assert(PyList_Check(py_centroids));
    assert(PyList_Check(py_data));
    assert(PyList_Check(py_indices));

    /* Convert Arrays from python to C */
    centroids_wrapper = py_list_to_array(py_centroids);
    data_points_wrapper = py_list_to_array(py_data);

    centroids = centroids_wrapper.pointers;
    data_points = data_points_wrapper.pointers;

    /* Main Algorithm - should call another func here which will call kmeans */
    result = kmeans(data_points, centroids, N, dim, K, max_iter);
    /* with the new interface:
    result = get_spk_clusters(data_points, centroids, N, dim, K, max_iter); */

    /* Print indices */
    for (i = 0; i < K; i++) {
        printf("%d", PyList_GET_ITEM(py_indices, i));
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
DataWrapper py_list_to_array(PyObject* py_list) {
    int n, m, i, j;
    double **result;
    double *p;
    PyObject *vector, *value;
    DataWrapper data_wrapper;

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
