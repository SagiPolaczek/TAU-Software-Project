/* Structures Declarations */
struct InputData {
    double **data;
    int dim;
    int data_count;
    double *container;
};
typedef struct InputData InputData;

struct CentroidsWrapper {
    double **centroids;
    double *container;
};
typedef struct CentroidsWrapper CentroidsWrapper;

struct Cluster{
    double *vector_sum;
    int *count;
};
typedef struct Cluster Cluster;

/* Functions Declarations */ 
InputData read_data();
CentroidsWrapper init_centroids(double **data_points, int K, int dim);
Cluster *init_clusters(int K, int dim);
int find_closest_centroid(double **centroids, double *data_point, int K, int d);
double compute_distance(double *u, double *v, int dim);
int update_centroids(double **centroids, Cluster *clusters, int K, int dim);
void add_datapoint_to_cluster(Cluster *clusters, int cluster_index, double *data_point, int dim);
void print_centroids(double **centroids, int K, int dim);
int isPositiveNumber(char* arg);