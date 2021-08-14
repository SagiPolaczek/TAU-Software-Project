#define PY_SSIZE_T_CLEAN
#include <Python.h>

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
