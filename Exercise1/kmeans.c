#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>



/* Structures Declerations */
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

/* Functions Decleration */ 
InputData read_data();
CentroidsWrapper init_centroids(double **data_points, int K, int dim);
Cluster *init_clusters(int K, int dim);
int find_closest_centroid(double **centroids, double *data_point, int K, int d);
double compute_distance(double *u, double *v, int dim);
int update_centroids(double **centroids, Cluster *clusters, int K, int dim);
void add_datapoint_to_cluster(Cluster *clusters, int cluster_index, double *data_point, int dim);
void print_centroids(double **centroids, int K, int dim);
int isPositiveNumber(char* arg);



/* Main Function */
int main(int argc, char *argv[]) 
{
    /* Variables Declerations */
    int K, max_iter;
    int data_count, dim = 0;
    int i;
    int seen_changes, count_iter, cluster_index;
    double **data_points, **centroids;
    double *data_point;
    double *data_container, *centroids_container;
    InputData input_data;
    CentroidsWrapper centroids_wrapper;
    Cluster *clusters;

    /* Manage Input */
    assert((argc == 2 || argc == 3) && "The program can have only 2 or 3 arguments!");
    
    assert(isPositiveNumber(argv[1]) && "Arguments must be a positive numbers!");
    K = atoi(argv[1]);
    assert(K > 0 && "K must be positive!");
    
    
    max_iter = 200;
    if (argc == 3) {
        assert(isPositiveNumber(argv[2]) && "Arguments must be a positive numbers!");
        max_iter = atoi(argv[2]);
        assert(max_iter >= 0 && "Max Iterations must be non-negative!");
    }

    /* Read data to a 2-Dimentional matrix () */
    input_data = read_data();

    data_points = input_data.data;
    dim = input_data.dim;
    data_count = input_data.data_count;
    data_container = input_data.container;
    
    /* The algorithm required K < N */ 
    assert(K < data_count && "K must be smaller than N!");

    /* Initalize centroids and clusters*/
    centroids_wrapper = init_centroids(data_points, K, dim);
    centroids = centroids_wrapper.centroids;
    centroids_container = centroids_wrapper.container;

    clusters = init_clusters(K, dim);

    /* Main Algorithm's Loop */
    count_iter = 0;
    seen_changes = 1;

    while ((seen_changes == 1) && (count_iter < max_iter)) {
        count_iter++;
        
        for (i = 0; i < data_count; i++) {
            data_point = data_points[i];
            cluster_index = find_closest_centroid(centroids, data_point, K, dim);
            add_datapoint_to_cluster(clusters, cluster_index, data_point, dim);            
        }

        seen_changes = update_centroids(centroids, clusters, K, dim);
    }

    /* Print centroids */
    print_centroids(centroids, K, dim);

    /* Free memory */
    free(data_points);
    free(data_container);
    free(centroids);
    free(centroids_container);
    free(clusters);

    return 0;
}

InputData read_data() {
    /* Variables Declerations */
    InputData data;
    double value;
    int data_count = 0, dim = 0;
    char c;
    double *p;
    double **data_points;
    int i, j;
    int first_round =1;

    /* Scan data from stdin TO COUNT dim and data_count */
    while (scanf("%lf%c", &value, &c) == 2)
    {
        /* Get dimention */
        if (first_round == 1) {
            dim++;
        }

        if (c == '\n') 
        {
            first_round = 0;
            /* Keep track on the vectors' amount */
            data_count++;
        }        
    }
    /* Rewind to the beginning of the data stream */
    rewind(stdin);

    /* Init a 2-dimentaional array for the data */ 
    p = calloc((data_count * dim), sizeof(double));
    assert(p != NULL);

    data_points = calloc(data_count, sizeof(double *));
    assert(data_points != NULL);

    for (i = 0; i < data_count; i++) {
        data_points[i] = p + i*dim;
    }
    
    /* Put the data from the stream into the Array */
    for (i = 0; i < data_count; i++) {
        for (j = 0; j < dim; j++) {
            scanf("%lf%c", &value, &c);
            data_points[i][j] = value;
        }
    }
    
    /* Load into the returning object */
    data.data_count = data_count;
    data.dim = dim;
    data.data = data_points;
    data.container = p;

    return data;
}

CentroidsWrapper init_centroids(double **data_points, int K, int dim){
    /* Variables Declerations */
    double *p;
    double **centroids;
    int i, j;
    CentroidsWrapper wrapper;

    /* Allocate space and define centroids */
    p = calloc((K * dim), sizeof(double));
    assert(p != NULL);

    centroids = calloc(K, sizeof(double *));
    assert(centroids != NULL);

    for (i = 0; i < K; i++) {
        centroids[i] = p + i*dim;
    }

    /* Put the first K data points into centroids */
    for (i = 0; i < K; i++) {
        for (j = 0; j < dim; j++) {
            centroids[i][j] = data_points[i][j];
        }
    }

    wrapper.centroids = centroids;
    wrapper.container = p;

    return wrapper;
}

Cluster *init_clusters(int K, int dim) {
    /* Variables Declerations */
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
    /* Variables Declerations */
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
    /* Variables Declerations */
    double distance = 0;
    int i = 0;

    /* Compute NORM^2 */
    for (i = 0; i < dim; i++) {
        distance += (u[i] - v[i]) * (u[i] - v[i]);
    }

    return distance;
}


int update_centroids(double **centroids, Cluster *clusters, int K, int dim) {  
    /* Variable Declerations */  
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
    /* Variables Declerations */
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

void print_centroids(double **centroids, int K, int dim) {
    /* Variables Declerations */
    double *centroid;
    double data_point;
    int i, j;

    /* Loop throughtout the centroids' data */
    for (i = 0; i < K; i++) {
        centroid = centroids[i];
        for (j = 0; j < dim; j++) {
            data_point = centroid[j];
            printf("%.4f", data_point);
            if (j < dim - 1) {
                /* seperate by comma adjecents data points */
                printf(",");
            }
        }
        
        /* seperate by new line adjecents centroids */
        printf("\n");
        
    }
}

int isPositiveNumber(char* arg) {
    int i, n;
    n = strlen(arg);
    for (i = 0; i < n; i++) {
        /* Check if arg[i] in {0,...,9} by comparing ASCII values */
        if (arg[i] < 48 || arg[i] > 57) {
            return 0;
        }
    }
    return 1;
}
