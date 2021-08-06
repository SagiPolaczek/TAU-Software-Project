#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


/* Structures Declarations */
struct DataWrapper {
    double **pointers;
    double *container;
};
typedef struct DataWrapper DataWrapper;

struct Cluster{
    double *vector_sum;
    int *count;
};
typedef struct Cluster Cluster;

/* Functions Declarations */
double** kmeans(double** data_points, double** centroids, int N, int dim, int K, int max_iter);
static PyObject* fit(PyObject* self, PyObject* args);
DataWrapper py_list_to_array(PyObject* py_list);
Cluster* init_clusters(int K, int dim);
int find_closest_centroid(double **centroids, double *data_point, int K, int d);
double compute_distance(double *u, double *v, int dim);
int update_centroids(double **centroids, Cluster *clusters, int K, int dim);
void add_datapoint_to_cluster(Cluster *clusters, int cluster_index, double *data_point, int dim);


/* Python <3 */
static PyObject* fit(PyObject* self, PyObject* args)
{
    /* Variables Declarations */
    PyObject *py_centroids, *py_data, *py_indices;
    PyObject *py_result, *py_vector;
    double **centroids, **data_points, **result;
    int N, dim, K, max_iter;
    int i, j;
    DataWrapper centroids_wrapper, data_points_wrapper;

    if (!PyArg_ParseTuple(args, "OOOiiii", &py_centroids, &py_data, &py_indices, &N, &dim, &K, &max_iter))
        return NULL;
    
    /* Validation */
    assert(PyList_Check(py_centroids));
    assert(PyList_Check(py_data));
    assert(PyList_Check(py_indices));

    /* Convert Arrays from python to C */
    centroids_wrapper = py_list_to_array(py_centroids);
    data_points_wrapper = py_list_to_array(py_data);

    centroids = centroids_wrapper.pointers;
    data_points = data_points_wrapper.pointers;

    /* Main Algorithm */
    result = kmeans(data_points, centroids, N, dim, K, max_iter);


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

    /* Free Memory */
    free(centroids_wrapper.container);
    free(data_points_wrapper.container);
    free(result);

    return py_result;
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



/* Main Function */
double** kmeans(double** data_points, double** centroids, int N, int dim, int K, int max_iter) 
{
    /* Variables Declarations */
    int i;
    int seen_changes, count_iter, cluster_index;
    double *data_point;
    Cluster *clusters;
        
    clusters = init_clusters(K, dim);

    /* Main Algorithm's Loop */
    count_iter = 0;
    seen_changes = 1;

    while ((seen_changes == 1) && (count_iter < max_iter)) {
        count_iter++;
        
        for (i = 0; i < N; i++) {
            data_point = data_points[i];
            cluster_index = find_closest_centroid(centroids, data_point, K, dim);
            add_datapoint_to_cluster(clusters, cluster_index, data_point, dim);            
        }

        seen_changes = update_centroids(centroids, clusters, K, dim);
    }

    /* Free memory */
    free(clusters);

    return centroids;
}

Cluster* init_clusters(int K, int dim) {
    /* Variables Declarations */
    Cluster *clusters;
    int i;
    int *count;
    double *array;

    /* Allocate space and Init clusters to default */
    clusters = calloc(K, sizeof(Cluster));
    assert(clusters != NULL);

    for (i = 0; i < K; i++) {
        /* Allocate array for each cluster */
        array = calloc(dim, sizeof(double));
        assert(array != NULL);
        
        /* Allocate pointer to count as a singleton */
        count = calloc(1, sizeof(int));
        assert(count != NULL);

        /* Attach */
        clusters[i].vector_sum = array;
        clusters[i].count = count;
    }

    return clusters;
}

int find_closest_centroid(double **centroids, double *data_point, int K, int dim){
    /* Variables Declarations */
    double min_distance;
    double curr_distance = 0;
    int i ,min_index = 0;
    
    min_distance = compute_distance(centroids[0], data_point, dim);
    /* Loop throughtout the cluster */
    for (i = 1; i < K; i++){
        curr_distance = compute_distance(centroids[i], data_point, dim);
        /* If we've found a closer centroid */
        if (curr_distance < min_distance) {
            min_distance = curr_distance;
            min_index = i;
        }
    }

    return min_index;
}


double compute_distance(double *u, double *v, int dim) {
    /* Variables Declarations */
    double distance = 0;
    int i = 0;

    /* Compute NORM^2 */
    for (i = 0; i < dim; i++) {
        distance += (u[i] - v[i]) * (u[i] - v[i]);
    }

    return distance;
}


int update_centroids(double **centroids, Cluster *clusters, int K, int dim) {  
    /* Variable Declarations */  
    Cluster cluster;
    double *cluster_vector;
    int cluster_count;
    double *centroid;
    double new_value;
    int i, j;
    int seen_changes = 0;

    /* Update centroids */
    for(i = 0 ; i < K ; i++) {
        cluster = clusters[i];
        cluster_vector = cluster.vector_sum;
        cluster_count = cluster.count[0];
        centroid = centroids[i];

        /* If cluster not empty */
        if(cluster_count > 0) {
            for(j = 0 ; j < dim ; j++) {
                new_value = (cluster_vector[j] / cluster_count);
                /* If there was a change */
                if(centroid[j] != new_value) {
                    centroid[j] = new_value;
                    seen_changes = 1;
                }
                /* Zerofy cluster's vector */
                cluster_vector[j] = 0;
            }
            /* Zerofy cluster's count */
            cluster.count[0] = 0;
        }
    }
    return seen_changes;
}

void add_datapoint_to_cluster(Cluster *clusters, int cluster_index,
                                         double *data_point, int dim){
    /* Variables Declarations */
    Cluster cluster;
    double *cluster_vector;
    int i;

    cluster = clusters[cluster_index];
    cluster_vector = cluster.vector_sum;


    /* Sum coordinate by coordinate */
    for (i = 0; i < dim; i++) {
        cluster_vector[i] += data_point[i];
    }
    
    /* Raise count by one */
    cluster.count[0] += 1;
}

/* Python Staff */
static PyMethodDef _methods[] = {
    {"fit", (PyCFunction)fit, METH_VARARGS, PyDoc_STR("Our SAVAGE Program")},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    return PyModule_Create(&_moduledef);
}
