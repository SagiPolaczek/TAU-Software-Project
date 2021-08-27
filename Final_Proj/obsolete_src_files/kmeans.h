#pragma once

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

typedef struct DataPoint {
    double *vector;
    int cluster_idx;
} DataPoint;


/* Functions Declarations */
void kmeans(double** data_points, int N, int dim, int K, int max_iter, double** centroids);
Cluster* init_clusters(int K, int dim);
int find_closest_centroid(double **centroids, double *data_point, int K, int d);
double compute_distance(double *u, double *v, int dim);
int update_centroids(double **centroids, Cluster *clusters, int K, int dim);
void add_datapoint_to_cluster(Cluster *clusters, int cluster_index, double *data_point, int dim);
void init_centroids(double **data_points, int K, int dim, double **centroids);
