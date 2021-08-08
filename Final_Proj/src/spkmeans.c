#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "spkmeans.h"
#include "debugger.h"

#define MEM_ALLOC_ERR "Fail to allocate memory."

/*
    TODO: 
    * Put the memory allocations into the functions themselves.
*/

/* CLI: k  goal  file_name */
int main(int argc, char *argv[]) {
    char *goal_string, *file_path;
    int k, N, i;
    goal goal;
    Graph graph = {0};
    double *ptr1;
    double *eigenvalues, **eigenvectors;

    LOG("\n----- DEBUGGING MODE ------\n\n");


    assert(argc == 4 && "Invalid arguments amount!\nOnly 3 arguments are allowed.");
    /* "You can assume that all input files and arguments are valid" - No need for further checks */
    
    k = atoi(argv[1]);
    goal_string = argv[2];
    goal = (int)(goal_string[0]);
    file_path = argv[3];
    
    /* Read the data into graph's {vertices, N, dim} */
    read_data(&graph, file_path);
    N = graph.n;
    

    compute_wam(&graph);

    /* If the goal is only to compute the WAM, we finished */
    if (goal == wam) {
        print_matrix(graph.weights, N, N);

        free_graph(&graph, wam);

        LOG("-- COMPLETE GOAL WAM --");
        return 42;
    }

    compute_ddg(&graph);
    
    if (goal == ddg) {
        print_vector_as_matrix(graph.degrees, N);

        free_graph(&graph, ddg);
        LOG("-- COMPLETE GOAL DDG --");
        return 42;
    }

    compute_lnorm(&graph);

    if (goal == lnorm) {
        print_matrix(graph.lnorm, N, N);

        free_graph(&graph, lnorm);
        LOG("-- COMPLETE GOAL LNORM --");
        return 42;
    }

    /* Allocate memory for the eigenvectors & eigenvalues */
    ptr1 = calloc((N * N), sizeof(double));
    assert(ptr1 != NULL && MEM_ALLOC_ERR);
    eigenvectors = calloc(N, sizeof(double*));
    assert(eigenvectors != NULL && MEM_ALLOC_ERR);
    for (i = 0; i < N; i++) {
        eigenvectors[i] = ptr1 + i*N;
    }

    eigenvalues = calloc(N, sizeof(double));
    assert(eigenvalues != NULL && MEM_ALLOC_ERR);

    compute_jacobi(&graph, eigenvectors, eigenvalues);

    if (goal == jacobi) {

        /* Print eigenvectors & eigenvalues */
        LOG("-- Eigenvectors:\n");
        print_matrix(eigenvectors, N, N);
        LOG("-- Eigenvalues:\n");
        print_matrix(&eigenvalues, 1, N);

        free(ptr1);
        free(eigenvectors);
        free(eigenvalues);

        free_graph(&graph, jacobi);

        LOG("-- COMPLETE GOAL JACOBI --");
        return 42;
    }

    /* goal == spk */

    /* Extra operations */
        
    /* Print result */

    free(ptr1);
    free(eigenvectors);
    free(eigenvalues);
    free_graph(&graph, spk);

    return 42;
}

/* Read the file and contain the data inplace in the graph */
void read_data(Graph *graph, char *file_path) {
    double value;
    char c;
    int first_round = 1;
    int data_count = 0, dim = 0;
    double *p;
    double **data_points;
    int i, j;
    FILE *ptr = fopen(file_path, "r");
    /* assert(ptr != NULL && "Could not load file"); */

    LOG("-- Reading Data --\n");

    /* Scan data from the file to *count* dim and data_count */
    while (fscanf(ptr, "%lf%c", &value, &c) == 2)
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

    rewind(ptr);

    /* Init a 2-dimentaional array for the data (vectors)*/ 
    p = calloc((data_count * dim), sizeof(double));
    assert(p != NULL && MEM_ALLOC_ERR);
    data_points = calloc(data_count, sizeof(double *));
    assert(data_points != NULL && MEM_ALLOC_ERR);

    for (i = 0; i < data_count; i++) {
        data_points[i] = p + i*dim;
    }

    /* Put the data from the stream into the Array */
    for (i = 0; i < data_count; i++) {
        for (j = 0; j < dim; j++) {
            fscanf(ptr, "%lf%c", &value, &c);
            data_points[i][j] = value;
        }
    }
    fclose(ptr);

    graph->vertices = data_points;
    graph->dim = dim;
    graph->n = data_count;    
    return;
}

void compute_wam(Graph *graph) {
    int i, j, N, dim;
    double weight, distance;
    double *v1, *v2;
    double **vertices = graph -> vertices;
    double **weights, *p;

    LOG("-- Computing WAM --\n");
    
    N = graph->n;
    dim = graph->dim;

    /* Allocate memory for the WAM */
    p = calloc((N * N), sizeof(double));
    assert(p != NULL && MEM_ALLOC_ERR);
    weights = calloc(N, sizeof(double*));
    assert(weights != NULL && MEM_ALLOC_ERR);
    for (i = 0; i < N; i++) {
        weights[i] = p + i*N;
    }

    graph->weights = weights;
    
    for (i = 0; i < N; i++) {
        v1 = vertices[i];
        for (j = 0; j < i; j++) {
            v2 = vertices[j];
            distance = compute_distance(v1, v2, dim);
            weight = exp((-1) * distance / 2);
            weights[i][j] = weight;
            weights[j][i] = weight;
        }
        weights[i][i] = 0;
    }
    LOG("-- COMPLETE COMPUTE WAM --\n");
    return;
}

double compute_distance(double *vec1, double *vec2, int dim) {
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
TO DO:
elaborate why we chose to represent the DDG with a vector instead of a matrix.
and how it affects the other functions
*/
void compute_ddg(Graph *graph) {
    int i;
    double deg;
    double **weights = graph -> weights;
    int N = graph -> n;
    double *ddg;

    LOG("-- Computing DDG --\n");

    /* Allocate memory for the DDG's "matrix" - represented by vector */
    ddg = calloc(N, sizeof(double*));
    assert(ddg != NULL && MEM_ALLOC_ERR);
    graph->degrees = ddg;

    for (i = 0; i < N; i++) {
        deg = compute_degree(weights, i, N);
        ddg[i] = deg;
    }
    return;
}

double compute_degree(double **weights, int v_idx, int n) {
    int i;
    double deg = 0;

    for (i = 0; i < n; i++) {
        deg += weights[v_idx][i];
    }

    return deg;
}

void compute_lnorm(Graph *graph) {
    double *invsqrt_d, *p;
    int i, j;
    double **weights = graph -> weights;
    double *degs = graph -> degrees;
    int N = graph -> n;
    double **lnorm;

    LOG("-- Computing LNORM --\n");

    /* Allocate memory for the Lnorm's matrix */
    p = calloc((N * N), sizeof(double));
    assert(p != NULL && MEM_ALLOC_ERR);
    lnorm = calloc(N, sizeof(double*));
    assert(lnorm != NULL && MEM_ALLOC_ERR);

    for (i = 0; i < N; i++) {
        lnorm[i] = p + i*N;
    }
    graph->lnorm = lnorm;

    invsqrt_d = inverse_sqrt_vec(degs, N);

    /* (invsqrt_d * W * invsqrt_d)
      = ((invsqrt_d * W) * invsqrt_d) */
    multi_vec_mat(invsqrt_d, weights, N, lnorm);
    multi_mat_vec(lnorm, invsqrt_d, N, lnorm);

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

double* inverse_sqrt_vec(double *vector, int n) {
    double *res;
    int i;

    /* Consider make it INPLACE */
    res = calloc(n, sizeof(double*));
    assert(res != NULL && MEM_ALLOC_ERR);
    
    for (i = 0; i < n; i++) {
        res[i] = 1 / (sqrt(vector[i]));
    }
    return res;
}

/* inplace */
void multi_vec_mat(double *vec, double **mat, int n, double **res) {
    int i, j;

    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            res[i][j] = vec[i] * mat[i][j];
        }
    }
}

void multi_mat_vec(double **mat, double *vec, int n, double **res) {
    int i, j;

    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            res[j][i] = vec[i] * mat[j][i];
        }
    } 
}

void compute_jacobi(Graph *graph, double **eign_vecs, double *eign_vals) {
    int is_not_diag = 1;
    int i, j, r, N;
    double max_entry, curr_entry;
    int max_row, max_col;
    double phi;
    double c, s;
    double **A_tag, **A;
    double *p;
    double s_sq, c_sq;
    double off_diff;
    const int MAX_ITER = 100;
    const double EPS = 0.001;
    int iter_count;

    LOG("-- Computing JACOBI --\n");

    N = graph->n;
    A = graph->lnorm;

    /* A_tag start as deep copy of A */
    p = calloc((N * N), sizeof(double));
    assert(p != NULL && MEM_ALLOC_ERR);
    A_tag = calloc(N, sizeof(double*));
    assert(A_tag != NULL && MEM_ALLOC_ERR);
    for (i = 0; i < N; i++) {
        A_tag[i] = p + i*N;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A_tag[i][j] = A[i][j];
        }
    }
    
    /* Initalize the Idendity Matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            eign_vecs[i][j] = (i == j);
        }
    }
    
    iter_count = 0;
    while(is_not_diag) {
        /* Pivot - Find the maximum absolute value A_ij */
        max_row = 0; max_col = 1; /* Assume n >= 2 */
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


        i = max_row;
        j = max_col;

        /* Obtain c,s (P) */
        phi = A[j][i];
        phi = 0.5*atan2(-2*A[i][j],A[j][j] - A[i][i]); /* Internet ??? */
        s = sin(phi);
        c = cos(phi);

        /* Relation between A and A_tag */
        i = max_row;
        j = max_col;
        for (r = 0; r < N; r++) {
            if ((r != i) && (r != j)) {
                A_tag[r][i] = c*A[r][i] - s*A[r][j];
                A_tag[r][j] = c*A[r][j] + s*A[r][i];
            }
        }
        c_sq = pow(c, 2);
        s_sq = pow(s, 2);
        A_tag[i][i] = c_sq*A[i][i] + s_sq*A[j][j] - 2*s*c*A[i][j];
        A_tag[j][j] = s_sq*A[i][i] + c_sq*A[j][j] + 2*s*c*A[i][j];
        A_tag[i][j] = (c_sq - s_sq)*A[i][j] + s*c*(A[i][i] - A[j][j]); /* => =0 */

        /* Update Eigenvectors' matrix */
        /* TO-DO: Check & fix! what if i=j? */ 
        for (r = 0; r < N; r++) {
                eign_vecs[r][j] = c*eign_vecs[r][i] - s*eign_vecs[r][j];
                eign_vecs[r][i] = c*eign_vecs[r][j] + s*eign_vecs[r][i];
        }

        /* Checking for Convergence */
        off_diff = 0;
        for (r = 0; r < N; r++) {
            if ((r != i) && (r != j)) {
                off_diff += pow(A[r][i], 2);
                off_diff += pow(A[r][j], 2);
                off_diff -= pow(A[r][i], 2);
                off_diff -= pow(A[r][j], 2);
            }
        }
        if (i != j) {
            off_diff += pow(A_tag[i][j], 2);
            off_diff -= pow(A[i][j], 2);
        }

        if (off_diff <= EPS || iter_count > MAX_ITER) {
            is_not_diag = 0;
            break;
        }
        
        /* 'Deep' update A = A' */
        for (r = 0; r < N; r++) {
            if ((r != i) && (r != j)) {
                A[r][i] = A_tag[r][i];
                A[r][j] = A_tag[r][j];
            }
        }
        A[i][i] = A_tag[i][i];
        A[j][j] = A_tag[j][j];
        A[i][j] = A_tag[i][j];

        iter_count += 1;
    }

    /* Update Eigenvalues */
    for (i = 0; i < N; i++) {
        eign_vals[i] = A_tag[i][i];
    }

    /* Free A_tag <3 */
    free(A_tag[0]);
    free(A_tag);
    
    return;
}

/* Sort the eigenvectors by the corresponding eigenvalues. INPLACE*/
void sort_by_eigen_values(double **vectors, double *values, int n) {
    /* Idea:
        qsort the eigenvalues. (fastest & easiest way (?)).
        Then order the eigenvectors to correspond the eigenvalues order.
        The last part will take O(n^2) but maybe should be pretty fast because n <= 1000.

        NOTE:
        We should change the matrix representation. 
        If so, we can achive O(1) time for swapping between two vectors instead of O(n).

        UPDATE (4.8 21:30)
        It might be better to use insertion sort! consider it.

        UPDATE (6.8 11:40)
        Another Idea!
        We can just sort the eigenvalues!
        If we sort only the eigenvalues, we can find the K easily.
        Then we'll copy only the K relevant eigenvectors by iterating over the matrix.
        The last operation bounded by O(n^2).

        UPDATE (6.8 13:34)
        Another Idea!
        We can use the qsort on the matrix itself and just define the comperator as we pleased!
        A comlexity of O(nlogn) without much effort! plus we get the fast implementation of the qsort.
        For this idea we shall use the transpose matrix, and make sure the corresponding values are
        swaps as their eigenvectors are being swap. 
    */

   double *values_copy;
   int i, j;
   double **vectors_T;
   const int NULL_VALUE = -42;
    
    /* Deep copy the eigenvalues */
    values_copy = calloc(n, sizeof(double*));
    assert(values_copy != NULL && MEM_ALLOC_ERR);
    for (i = 0; i < n; i++) {
        values_copy[i] = values[i];
    }

    /* Sort the values */   
    qsort(values, n, sizeof(int), cmpfunc);

    /* Copy the vectors after transpose into an array which allow us to swap rows (columns) in O(1) */
    vectors_T = calloc(n, sizeof(int *));
    assert(vectors_T != NULL && MEM_ALLOC_ERR);
    for (i = 0; i < n; i++) {
        vectors_T[i] = calloc(n, sizeof(int));
        assert(vectors_T[i] != NULL && MEM_ALLOC_ERR);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            vectors_T[j][i] = vectors[i][j];
        }
    }


    /* Order the rows (columns) of vectors_T to correspond the *sorted* values */
    for (i = 0; i < n; i++) {
        int curr_sorted_value = values[i];
        for (j = 0; j < n; j++) {
            int curr_value = values_copy[j];
            if (curr_value == curr_sorted_value) {
                /* need to swap */
                /* TODO: consider extract to outer function */
                double *tmp = vectors_T[i];
                vectors_T[i] = vectors_T[j];
                vectors_T[j] = tmp;

                /* mark the value as visited */
                values_copy[j] = NULL_VALUE;
                /* moved to the next sorted value */
                continue;
            }
        }
    }

    /* copy the sorted rows as sorted columns */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            vectors[j][i] = vectors_T[i][j];
        }
    }
    
    for (i = 0; i < n; i++) {
        free(vectors_T[i]);
    }
    free(vectors_T);

}

/* Define a comparator for the qsort function.
   Credit: 'tutorialspoint' */
int cmpfunc (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void free_graph(Graph *graph, goal goal) {
    LOG("-- FREE GRAPH --\n");
    free((graph->vertices)[0]);
    free(graph->vertices);
    free((graph->weights)[0]);
    free(graph->weights);

    if (goal == wam) {return;}

    free(graph->degrees);

    if (goal == ddg) {return;}

    free((graph->lnorm)[0]);
    free(graph->lnorm);

    return;
}

void print_matrix(double **mat, int rows, int cols) {
    double *vec;
    double val;
    int i, j;

    for (i = 0; i < rows; i++) {
        vec = mat[i];
        for (j = 0; j < cols; j++) {
            val = vec[j];
            printf("%.4f", val);
            if (j < cols - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

/* Final function for operate the SPK.
   We assume that jacobi was computed (and all it's dependencies) */

void finish_spk() {
    
    /*  If k doesn't provided (k==0):
            Determine k */

    /*  Obtain the first k eigenvectors u1, ..., uk
        Form the matrix U(nxk) which u_i is the i'th column */

    /*  Form T from U by renormalizing each row to have the unit length */

    /*  Traet each row of T as a point in R^k.
        cluster them into k clusters with the k-means */

    /*  Assign the original points to the clusters */
    
}
