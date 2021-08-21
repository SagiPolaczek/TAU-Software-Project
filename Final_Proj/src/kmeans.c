#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "kmeans.h"

/* Main Function */
extern double** kmeans(double** data_points, double** centroids, int N, int dim, int K, int max_iter) 
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
