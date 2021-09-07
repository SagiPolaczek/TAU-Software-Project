#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "spkmeans.h"

#define ERR_MSG "An Error Has Occured"

const int MAX_ITER_KMEANS = 300;
const int MAX_ITER_JACOBI = 100;

/*
    Command Line Interface in the following format:
    argv = [K, goal, file_name], argc = 4
*/
int main(int argc, char *argv[])
{
    char *goal_string, *file_path;
    double **data_points, **centroids;
    int K, N;
    goal goal;
    Graph graph = {0};

    if (argc != 4)
    {
        printf("Only 3 Arguments are allowed.\n");
        exit(1);
    }

    K = atoi(argv[1]);
    goal_string = argv[2];
    goal = (int)(goal_string[0]);
    file_path = argv[3];

    /* Read the data into graph's {vertices, N, dim} */
    read_data(&graph, file_path);
    N = graph.N;

    /* goal == wam / ddg / lnorm / jacobi */
    if (goal != spk)
    {
        compute_by_goal(&graph, goal);
    }
    /* goal == spk */
    else
    {
        data_points = init_spk_datapoints(&graph, &K);
        print_matrix(data_points, N, K);
        centroids = calloc_matrix(K, K);
        init_centroids(data_points, K, K, centroids);
        print_matrix(centroids, K, K);
        kmeans(data_points, N, K, K, MAX_ITER_KMEANS, centroids);

        print_matrix(centroids, K, K);

        free_graph(&graph, goal);
        free_matrix(centroids);
        free_matrix(data_points);
    }

    return 42;
}

/*
    Read the file and store the data in the graph
*/
void read_data(Graph *graph, char *file_path)
{
    double value;
    char c;
    int first_round = 1;
    int data_count = 1, dim = 0;
    double **data_points;
    int i, j;
    FILE *ptr = fopen(file_path, "r");
    if (ptr == NULL) {
        printf(ERR_MSG);
        exit(1);
    }

    /* Scan data from the file to *count* dim and data_count */
    while (fscanf(ptr, "%lf%c", &value, &c) == 2)
    {
        /* Get dimention */
        if (first_round == 1) {
            dim++;
        }

        if (c == '\n') {
            first_round = 0;
            /* Keep track on the vectors' amount */
            data_count++;
        }
    }

    rewind(ptr);

    /* Init a 2-dimentaional array for the data (vectors)*/
    data_points = calloc_matrix(data_count, dim);

    /* Put the data from the stream into the Array */
    for (i = 0; i < data_count; i++)
    {
        for (j = 0; j < dim; j++)
        {
            fscanf(ptr, "%lf%c", &value, &c);
            data_points[i][j] = value;
        }
    }
    fclose(ptr);

    graph->vertices = data_points;
    graph->dim = dim;
    graph->N = data_count;
    return;
}

/*
    Compute the Weighted Adjacency Matrix and store it in the graph.
*/
void compute_wam(Graph *graph)
{
    int i, j, N, dim;
    double weight, distance;
    double *v1, *v2;
    double **vertices = graph->vertices;
    double **weights;

    N = graph->N;
    dim = graph->dim;

    /* Allocate memory for the WAM */
    weights = calloc_matrix(N, N);
    graph->weights = weights;

    for (i = 0; i < N; i++)
    {
        v1 = vertices[i];

        for (j = 0; j < i; j++)
        {
            v2 = vertices[j];
            distance = compute_distance_spk(v1, v2, dim);
            weight = exp((-1) * (distance / 2));
            weights[i][j] = weight;
            weights[j][i] = weight;
        }
        weights[i][i] = 0;
    }
    return;
}

/*
    Compute the distance between the two vector.
    The vectors have 'dim' as dimention.
 */
double compute_distance_spk(double *vec1, double *vec2, int dim)
{
    /* Variables Declarations */
    int i;
    double distance = 0, res;

    /* Compute NORM */
    for (i = 0; i < dim; i++) {
        distance += pow((vec1[i] - vec2[i]), 2);
    }

    res = sqrt(distance);
    return res;
}

/*
    Compute the Diagonal Degree Matrix and store it in the graph.

    NOTE:
    We chose to store the matrix in an array data structure because
    all the entries are zero except from the diagonal. So a N dimentional
    vector is sufficient.
    It required us to define special functions for multiplication and printing,
    where we refer DDG's vector as a matrix.
*/
void compute_ddg(Graph *graph)
{
    int i, N;
    double deg;
    double **weights, *ddg;

    weights = graph->weights;
    N = graph->N;

    /* Allocate memory for the DDG's "matrix" - represented by vector */
    ddg = calloc_vector(N);
    graph->degrees = ddg;

    for (i = 0; i < N; i++) {
        deg = compute_degree(weights, i, N);
        ddg[i] = deg;
    }
    return;
}

/*
    Compute the degree of a vertex in the graph.
*/
double compute_degree(double **weights, int v_idx, int n)
{
    int i;
    double deg = 0;

    for (i = 0; i < n; i++) {
        deg += weights[v_idx][i];
    }

    return deg;
}

/*
    Compute the Normalized Graph Laplacian (Lnorm matrix).
*/
void compute_lnorm(Graph *graph)
{
    double *invsqrt_d;
    int i, j;
    double **weights = graph->weights;
    double *degs = graph->degrees;
    int N = graph->N;
    double **lnorm;

    /* Allocate memory for the Lnorm's matrix */
    lnorm = calloc_matrix(N, N);
    graph->lnorm = lnorm;

    invsqrt_d = calloc_vector(N);

    inverse_sqrt_vec(degs, N, invsqrt_d);

    multi_vec_mat_vec(invsqrt_d, weights, N, lnorm);

    free(invsqrt_d);

    /* res = I - res */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            lnorm[i][j] = (-1) * lnorm[i][j];
        }
        lnorm[i][i] += 1;
    }
    return;
}

/* Compute inv_sqrt_vec when:
   inv_sqrt_vec = (vector)^(-0.5) */
void inverse_sqrt_vec(double *vector, int N, double *inv_sqrt_vec)
{
    int i;

    for (i = 0; i < N; i++) {
        inv_sqrt_vec[i] = pow(vector[i], -0.5);
    }
    return;
}

/*
    Compute multiplication of: ' VECTOR * MATRIX * VECTOR'
    when the 'VECTOR' is a diagonal matrix A such that A[i][i]=V[i]
*/
void multi_vec_mat_vec(double *vec, double **mat, int n, double **res)
{
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i][j] = vec[i] * mat[i][j] * vec[j];
        }
    }
}

/*
    Compute Jacobi as has been described in the task.
    When: A_tag := A', eigenvectors := V (Corresponding to the task).
*/
void compute_jacobi(double **A, int N, double **eigenvectors, double *eigenvalues)
{
    int is_not_diag = 1;
    int i, j, r;
    int max_row, max_col;
    int sign_theta;
    double theta, c, s, t;
    double max_entry, curr_entry;
    double **A_tag, **eigenvectors_deepcopy;
    double off_diff;
    int iter_count;
    const double EPS = pow(10, -15);

    /* A_tag start as a deep copy of A */
    A_tag = calloc_matrix(N, N);
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            A_tag[i][j] = A[i][j];
        }
    }

    eigenvectors_deepcopy = calloc_matrix(N, N);
    init_idendity_matrix(N, eigenvectors);
    init_idendity_matrix(N, eigenvectors_deepcopy);

    iter_count = 1;
    while (is_not_diag) {
        /* Pivot - Find the maximum absolute value A_ij */
        max_row = 0;
        max_col = 1; /* Assume n >= 2 */
        max_entry = fabs(A[max_row][max_col]);

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                curr_entry = fabs(A[i][j]);
                if ((i != j) && (curr_entry > max_entry)) {
                    max_entry = curr_entry;
                    max_row = i;
                    max_col = j;
                }
            }
        }
        
        /* Update variables to match the task's description */
        i = max_row;
        j = max_col;

        /* Obtain c,s (P matrix) */
        theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
        sign_theta = get_sign(theta);
        t = sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / sqrt(pow(t, 2) + 1);
        s = t * c;

        update_A_tag(A, A_tag, N, i, j, c, s);
        update_eigenvectors(eigenvectors, eigenvectors_deepcopy, N, i, j, c, s);

        /* Compute the difference: off(A)^2 - off(A')^2 */
        off_diff = compute_off_diagonal_difference(A, A_tag, N, i, j);

        update_A(A, A_tag, N, i, j);

        if (off_diff <= EPS || iter_count >= MAX_ITER_JACOBI) {
            is_not_diag = 0; /* Is diagonal */
            break;
        }

        iter_count += 1;
    } /* While is not diagonal */

    /* Update Eigenvalues */
    for (r = 0; r < N; r++) {
        eigenvalues[r] = A[r][r];
    }

    /* Free <3 */
    free_matrix(A_tag);
    free_matrix(eigenvectors_deepcopy);

    return;
}

void update_eigenvectors(double **eigenvectors, double **eigenvectors_deepcopy,
                        int N, int i, int j, double c, double s)
{
    int r;
    /* Update Eigenvectors' matrix */
    for (r = 0; r < N; r++) {
        eigenvectors[r][i] = c * eigenvectors_deepcopy[r][i] - s * eigenvectors_deepcopy[r][j];
        eigenvectors[r][j] = c * eigenvectors_deepcopy[r][j] + s * eigenvectors_deepcopy[r][i];
    }

    /* Update the deep copy of the eigenvectors*/
    for (r = 0; r < N; r++) {
        eigenvectors_deepcopy[r][i] = eigenvectors[r][i];
        eigenvectors_deepcopy[r][j] = eigenvectors[r][j];
    }
}

void update_A_tag(double **A, double **A_tag, int N, int i, int j, double c, double s)
{
    double c_squared, s_squared;
    int r;
    for (r = 0; r < N; r++) {
        if ((r != i) && (r != j)) {
            A_tag[r][i] = c * A[r][i] - s * A[r][j];
            A_tag[r][j] = c * A[r][j] + s * A[r][i];
            A_tag[i][r] = A_tag[r][i];
            A_tag[j][r] = A_tag[r][j];
        }
    }
    c_squared = pow(c, 2);
    s_squared = pow(s, 2);
    A_tag[i][i] = c_squared * A[i][i] + s_squared * A[j][j] - 2 * s * c * A[i][j];
    A_tag[j][j] = s_squared * A[i][i] + c_squared * A[j][j] + 2 * s * c * A[i][j];
    A_tag[i][j] = 0; /* (c_squared - s_squared)*A[i][j] + s*c*(A[i][i] - A[j][j]) => =0 */
    A_tag[j][i] = 0;
}

void update_A(double **A, double **A_tag, int N, int i, int j)
{
    int r;
    for (r = 0; r < N; r++){
        A[r][i] = A_tag[r][i];
        A[r][j] = A_tag[r][j];
        A[i][r] = A_tag[r][i];
        A[j][r] = A_tag[r][j];
    }
}

double compute_off_diagonal_difference(double **A, double **A_tag, int N, int i, int j)
{
    double off_diff = 0;
    int r;

    for (r = 0; r < N; r++){
        off_diff += pow(A[r][i], 2);
        off_diff += pow(A[r][j], 2);
        off_diff -= pow(A_tag[r][i], 2);
        off_diff -= pow(A_tag[r][j], 2);
        }

    off_diff -= pow(A[i][i], 2);
    off_diff -= pow(A[j][j], 2);
    off_diff += pow(A_tag[i][i], 2);
    off_diff += pow(A_tag[j][j], 2);

    return off_diff;
}

/* 
    Initialize 'matrix' as the idendity matrix I 
*/
void init_idendity_matrix(int N, double** matrix)
{
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
        {
            matrix[i][j] = (i == j);
        }
    }
}


/*
    Define a comperator for the qsort function.
    SOURCE: 'https://stackoverflow.com/questions/20584499/why-qsort-from-stdlib-doesnt-work-with-double-values-c'
*/
int cmp_func(const void *a, const void *b)
{
    if (*(double *)a > *(double *)b)
        return 1;
    else if (*(double *)a < *(double *)b)
        return -1;
    else
        return 0;
}

/*
    Print the matrix with the given dimention in the right format.
*/
void print_matrix(double **mat, int rows, int cols)
{
    double *vec;
    double val;
    int i, j;

    for (i = 0; i < rows; i++) {
        vec = mat[i];
        for (j = 0; j < cols; j++) {
            val = vec[j];
            /* Avoiding '-0.0000' */
            if (val < 0 && val > -0.00005) {
                val = 0;
            }
            printf("%.4f", val);
            if (j < cols - 1) {
                printf(",");
            }
        }
        if (i < rows - 1) {
            printf("\n");
        }
    }
}

/*
    Print the matrix with the given dimention in the right format,
    when PRINTED[i][j] = MATRIX[j][i].

    NOTE:
    We define a special function for that because we didn't want to spend time on rotating the matrix.
    Its a tradeoff of: Time Vs Additional code and elegant. We chose time :) .
*/
void print_transpose_matrix(double **mat, int rows, int cols)
{
    double val;
    int i, j;

    for (i = 0; i < cols; i++) {
        for (j = 0; j < rows; j++) {
            val = mat[j][i];
            if (val < 0 && val > -0.00005) {
                val = 0;
            }
            printf("%.4f", val);
            if (j < rows - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

/*
    Print a vector as a digonal matrix.
    (used as we elaborated in the 'compute_ddg' function)
*/
void print_vector_as_matrix(double *diag, int n)
{
    double val;
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                val = diag[i];
            }
            else
            {
                val = 0;
            }
            printf("%.4f", val);
            if (j < n - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

/*
    Preforms steps 1-5 of the algorithm.
    
    TODO:
    elaborate.
*/
double **init_spk_datapoints(Graph *graph, int *K)
{
    double **eigenvectors, *eigenvalues, *eigenvalues_sorted;
    int N, i;
    double **U, **T;

    compute_wam(graph);
    compute_ddg(graph);
    compute_lnorm(graph);

    /* Allocate memory for the eigenvectors & eigenvalues */
    N = graph->N;
    eigenvectors = calloc_matrix(N, N);
    eigenvalues = calloc_vector(N);

    compute_jacobi(graph->lnorm, N, eigenvectors, eigenvalues);

    /* Deep copy the eigenvalues */
    eigenvalues_sorted = calloc_vector(N);

    for (i = 0; i < N; i++) {
        eigenvalues_sorted[i] = eigenvalues[i];
    }

    /* Sort eigenvalues_sorted */
    qsort(eigenvalues_sorted, N, sizeof(double), cmp_func);

    /*  If K doesn't provided (k==0) determine k */
    if (*K == 0) {
        *K = get_heuristic(eigenvalues_sorted, N);
    }

    /*  Obtain the first (ordered!) k eigenvectors u1, ..., uk
        Form the matrix U(nxk) which u_i is the i'th column */

    /* Allocate memory for the U matrix. nxk. */
    U = calloc_matrix(N, (*K));

    form_U(U, eigenvectors, eigenvalues, eigenvalues_sorted, N, *K);

    /*  Form T from U by renormalizing each row to have the unit length */
    T = U;
    form_T(T, N, *K);

    free(eigenvalues_sorted);
    free(eigenvalues);
    free_matrix(eigenvectors);

    return T;
}

/*
    TODO:
    elaborate.
*/
int get_heuristic(double *eigenvalues, int N)
{
    double max_delta = 0;
    int max_idx = ((N / 2) + 1 < N - 2) ? (N / 2) + 1 : N - 2; /* i is bounded as shown in the pdf */
    int i, K = 0;
    double val1, val2, curr_delta;

    for (i = 0; i < max_idx; i++) {
        val1 = eigenvalues[i];
        val2 = eigenvalues[i + 1];
        curr_delta = val2 - val1;
        if (curr_delta > max_delta) {
            max_delta = curr_delta;
            K = i + 1;
        }
    }
    return K;
}

/*
    Form U, the matrix containing the first K eigenvectors.
    NOTE:
    By "the first K eigenvectors" we refer to the eigenvectors that correspond 
    to the K smallest eigenvalues.
*/
void form_U(double **U, double **eigenvectors, double *eigenvalues, double *eigenvalues_sorted, int N, int K)
{
    double curr_val, curr_sorted_val;
    const int NULL_VAL = -42; /* define null value when known that eigenvalues are positive */
    int i, j, r;

    for (i = 0; i < K; i++) {
        curr_sorted_val = eigenvalues_sorted[i];
        for (j = 0; j < N; j++) {
            curr_val = eigenvalues[j];
            if (curr_val == curr_sorted_val)
            {
                /* place the j'th vector in the U's i'th column */
                for (r = 0; r < N; r++) {
                    U[r][i] = eigenvectors[r][j]; /* assume the eigenvectors arrage as columns - DOUBLE CHECK IT */
                }

                /* mark the eigenvalue as "visited" */
                eigenvalues[j] = NULL_VAL;

                /* move to the next sorted eigenvalue */
                break;
            }
        }
    }
}

/*
    Form the matrix T from U by renormalizing each of U's rows.
    NOTE:
    T,U naming is corresponding to the task's description.
*/
void form_T(double **U, int N, int K)
{
    int i, j;
    double vec_len;
    double *vec, *zero_vec;

    /* init a vector of zeros */
    zero_vec = calloc_vector(K);

    for (i = 0; i < N; i++) {
        vec = U[i];
        /* compute the row vector's length */
        vec_len = compute_distance_spk(vec, zero_vec, K);

        for (j = 0; j < K; j++) {
            vec[j] = (vec[j] / vec_len);
        }
    }

    free(zero_vec);
}

/*
    Allocating memory for a vector.
*/
double *calloc_vector(int size)
{
    double *array;
    array = calloc(size, sizeof(double));
    my_assert(array != NULL);

    return array;
}

/*
    Allocating memory for a matrix.
*/
double **calloc_matrix(int rows, int cols)
{
    double *ptr, **array;
    int i;

    ptr = calloc((rows * cols), sizeof(double));
    array = calloc(rows, sizeof(double *));
    my_assert(array != NULL || ptr != NULL);
    for (i = 0; i < rows; i++)
    {
        array[i] = ptr + i * cols;
    }

    return array;
}

/*
    Free a memory of a matrix.
*/
void free_matrix(double **matrix)
{
    free(matrix[0]);
    free(matrix);
}

/*
    Free the memory of the graph.
*/
void free_graph(Graph *graph, goal goal)
{

    free_matrix(graph->vertices);
    if (goal == jacobi) {
        return;
    }

    free_matrix(graph->weights);
    if (goal == wam) {
        return;
    }

    free(graph->degrees);

    if (goal == ddg) {
        return;
    }

    free_matrix(graph->lnorm);
}

/*
    Return the sign of a double.
*/
int get_sign(double d)
{
    if (d >= 0)
        return 1;
    else
        return -1;
}

/*
    Custom assert fucntion to satisfy the assignment's requirements.
*/
void my_assert(int status)
{
    if (status == 0) {
        printf(ERR_MSG);
        exit(1);
    }
}

void compute_by_goal(Graph *graph, goal goal)
{
    double *eigenvalues, **eigenvectors, **A;
    int N = graph->N;

    if (goal == jacobi) {

        /* Allocate memory for the eigenvectors & eigenvalues */
        eigenvectors = calloc_matrix(N, N);
        eigenvalues = calloc_vector(N);

        A = graph->vertices;
        compute_jacobi(A, N, eigenvectors, eigenvalues);

        /* Print eigenvectors & eigenvalues */
        print_matrix(&eigenvalues, 1, N);

        /* Transpose to print the eigenvectors as rows */
        print_transpose_matrix(eigenvectors, N, N);

        free_matrix(eigenvectors);
        free(eigenvalues);

        free_graph(graph, jacobi);

        return;
    }

    compute_wam(graph);

    /* If the goal is only to compute the WAM, we finished */
    if (goal == wam) {
        print_matrix(graph->weights, N, N);

        free_graph(graph, wam);

        return;
    }

    compute_ddg(graph);

    if (goal == ddg) {
        print_vector_as_matrix(graph->degrees, N);

        free_graph(graph, ddg);
        return;
    }

    compute_lnorm(graph);

    if (goal == lnorm) {
        print_matrix(graph->lnorm, N, N);

        free_graph(graph, lnorm);
        return;
    }
}

/* Main Function */
void kmeans(double **data_points, int N, int dim, int K, int max_iter, double **centroids)
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

    while ((seen_changes == 1) && (count_iter < max_iter))
    {
        count_iter++;

        for (i = 0; i < N; i++)
        {
            data_point = data_points[i];
            cluster_index = find_closest_centroid(centroids, data_point, K, dim);
            add_datapoint_to_cluster(clusters, cluster_index, data_point, dim);
        }

        seen_changes = update_centroids(centroids, clusters, K, dim);
    }

    /* Free memory */
    free(clusters);

    return;
}

Cluster *init_clusters(int K, int dim)
{
    /* Variables Declarations */
    Cluster *clusters;
    int i;
    int *count;
    double *array;

    /* Allocate space and Init clusters to default */
    clusters = calloc(K, sizeof(Cluster));
    assert(clusters != NULL);

    for (i = 0; i < K; i++)
    {
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

int find_closest_centroid(double **centroids, double *data_point, int K, int dim)
{
    /* Variables Declarations */
    double min_distance;
    double curr_distance = 0;
    int i, min_index = 0;

    min_distance = compute_distance(centroids[0], data_point, dim);
    /* Loop throughtout the cluster */
    for (i = 1; i < K; i++)
    {
        curr_distance = compute_distance(centroids[i], data_point, dim);
        /* If we've found a closer centroid */
        if (curr_distance < min_distance)
        {
            min_distance = curr_distance;
            min_index = i;
        }
    }

    return min_index;
}

double compute_distance(double *u, double *v, int dim)
{
    /* Variables Declarations */
    double distance = 0;
    int i = 0;

    /* Compute NORM^2 */
    for (i = 0; i < dim; i++)
    {
        distance += (u[i] - v[i]) * (u[i] - v[i]);
    }

    return distance;
}

int update_centroids(double **centroids, Cluster *clusters, int K, int dim)
{
    /* Variable Declarations */
    Cluster cluster;
    double *cluster_vector;
    int cluster_count;
    double *centroid;
    double new_value;
    int i, j;
    int seen_changes = 0;

    /* Update centroids */
    for (i = 0; i < K; i++)
    {
        cluster = clusters[i];
        cluster_vector = cluster.vector_sum;
        cluster_count = cluster.count[0];
        centroid = centroids[i];

        /* If cluster not empty */
        if (cluster_count > 0)
        {
            for (j = 0; j < dim; j++)
            {
                new_value = (cluster_vector[j] / cluster_count);
                /* If there was a change */
                if (centroid[j] != new_value)
                {
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
                              double *data_point, int dim)
{
    /* Variables Declarations */
    Cluster cluster;
    double *cluster_vector;
    int i;

    cluster = clusters[cluster_index];
    cluster_vector = cluster.vector_sum;

    /* Sum coordinate by coordinate */
    for (i = 0; i < dim; i++)
    {
        cluster_vector[i] += data_point[i];
    }

    /* Raise count by one */
    cluster.count[0] += 1;
}

void init_centroids(double **data_points, int K, int dim, double **centroids)
{
    /* Variables Declarations */
    int i, j;

    /* Put the first K data points into centroids */
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < dim; j++)
        {
            centroids[i][j] = data_points[i][j];
        }
    }

    return;
}
