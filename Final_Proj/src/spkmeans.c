#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "spkmeans.h"
#include "debugger.h"

#define _MEM_ALLOC_ERR "Fail to allocate memory."




typedef enum goal { spk = (int)'s', 
                    wam = (int)'w',
                    ddg = (int)'d', 
                    lnorm = (int)'l', 
                    jacobi = (int)'j' } goal;

/* CLI: k  goal  file_name */
int main(int argc, char *argv[]) {
    char *goal_string, *file_path;
    int k;
    goal goal;
    /* Graph graph = {0}; */


    LOG("\n----- DEBUGGING MODE ------\n");


    assert(argc == 4 && "Invalid arguments amount!\nOnly 3 arguments are allowed.");
    /* "You can assume that all input files and arguments are valid" - No need for further checks */
    
    k = atoi(argv[1]);

    goal_string = argv[2];
    goal = (int)(goal_string[0]);

    file_path = argv[3];

    

    /* read_data(file_path, &graph); */ 
    


    switch (goal)
    {
        case spk:
        {

        } break;
        
        case wam:
        {

        } break;

        case ddg:
        {

        } break;

        case lnorm:
        {

        } break;

        case jacobi:
        {

        } break;
    }

    printf("%d%s", argc, argv[0]); /* For compilation purposes only. Remove when finshed */
    return 0;
}

/* Read the file and contain the data inplace in the graph */
void read_data(char *file_path, Graph *graph) {
    double value;
    char c;
    int first_round = 1;
    int data_count = 0, dim = 0;
    double *p;
    double **data_points;
    int i, j;


    FILE *ptr = fopen(file_path, "r");
    assert(ptr != NULL && "Could not load file");

    /* Scan data from the file to *count* dim and data_count */
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
    
    rewind(ptr);

    /* Init a 2-dimentaional array for the data (Vectors)*/ 
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

    graph->vertices = data_points;
    graph->dim = dim;
    graph->n = data_count;
}

double **compute_wam(Graph *g) {
   /* int i, j, n;
    double weight, distance;
    Vector *v1, *v2;
    Vector **vertices = g -> vertices;
    double **vertices = g -> vertices;

    */

    double **weights = g -> weights;

   /*  n = g -> n;
    for (i = 0; i < n; i++) {
        v1 = vertices[i];
        for (j = 0; j < i; j++) {
            v2 = vertices[j];
            distance = compute_distance(v1, v2);
            weight = exp((-1) * distance / 2);
            weights[i][j] = weight;
            weights[j][i] = weight;
        }
        weights[i][i] = 0;
    } */ 
    return weights;
}

double compute_distance(Vector *vec1, Vector *vec2) {
    /* Variables Declarations */
    int i;
    double distance = 0, res;
    int dim = vec1 -> dim;
    
    double *u = vec1 -> values;
    double *v = vec2 -> values;

    /* Compute NORM */
    for (i = 0; i < dim; i++) {
        distance += pow((u[i] - v[i]), 2);
    }

    res = sqrt(distance);
    return res;
}

/* 
TO DO:
elaborate why we chose to represent the DDG with a vector instead of a matrix.
and how it affects the other functions
*/
double *compute_ddg(Graph *g) {
    int i;
    double deg;
    double **weights = g -> weights;
    int n = g -> n;
    double *degs = g -> degrees;

    for (i = 0; i < n; i++) {
        deg = compute_degree(weights, i, n);
        degs[i] = deg;
    }
    return degs;
}

double compute_degree(double **weights, int v_idx, int n) {
    int i;
    double deg = 0;

    for (i = 0; i < n; i++) {
        deg += weights[v_idx][i];
    }

    return deg;
}

double** compute_lnorm(Graph *g) {
    double *invsqrt_d, **res, *p;
    int i, j;
    double **weights = g -> weights;
    double *degs = g -> degrees;
    int n = g -> n;

    p = calloc((n * n), sizeof(double));
    assert(p != NULL && _MEM_ALLOC_ERR);
    res = calloc(n, sizeof(double*));
    assert(res != NULL && _MEM_ALLOC_ERR);

    for (i = 0; i < n; i++) {
        res[i] = p + i*n;
    }

    invsqrt_d = inverse_sqrt_vec(degs, n);

    /* (invsqrt_d * W * invsqrt_d)
      = ((invsqrt_d * W) * invsqrt_d) */
    multi_vec_mat(invsqrt_d, weights, n, res);
    multi_mat_vec(res, invsqrt_d, n, res);

    /* res = I - res */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i][j] = (-1) * res[i][j];
        }
        res[i][i] += 1;
    }
    return res;
}

double* inverse_sqrt_vec(double *vector, int n) {
    double *res;
    int i;

    /* Consider make it INPLACE */
    res = calloc(n, sizeof(double*));
    assert(res != NULL && _MEM_ALLOC_ERR);
    
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

void jacobi_alg(double **A, int n, double **eign_vecs, double *eign_vals) {
    int is_not_diag = 1;
    int i, j, r;
    double max_entry, curr_entry;
    int max_row, max_col;
    double phi;
    double c, s;
    double **A_tag;
    double *p;
    double s_sq, c_sq;
    double off_diff;
    const int MAX_ITER = 100;
    const double EPS = 0.001;
    int iter_count;


    /* A_tag start as deep copy of A */
    p = calloc((n * n), sizeof(double));
    assert(p != NULL && _MEM_ALLOC_ERR);
    A_tag = calloc(n, sizeof(double*));
    assert(A_tag != NULL && _MEM_ALLOC_ERR);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A_tag[i][j] = A[i][j];
        }
    }

    /* Initalize the Idendity Matrix */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            eign_vecs[i][j] = (i == j);
        }
    }
    
    iter_count = 0;
    while(is_not_diag) {
        /* Pivot - Find the maximum absolute value A_ij */
        max_row = 0; max_col = 1; /* Assume n >= 2 */
        max_entry = fabs(A[max_row][max_col]);

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                curr_entry = fabs(A[i][j]);
                if ((i != j) && (curr_entry > max_entry)) {
                    max_entry = curr_entry;
                    max_row = i;
                    max_col = j;
                }
            }
        }

        /* Obtain c,s (P) */
        phi = 0.5*atan2(-2*A[i][j],A[j][j] - A[i][i]); /* Internet ??? */
        s = sin(phi);
        c = cos(phi);

        /* Relation between A and A_tag */
        i = max_row;
        j = max_col;
        for (r = 0; r < n; r++) {
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
        for (r = 0; r < n; r++) {
                eign_vecs[r][j] = c*eign_vecs[r][i] - s*eign_vecs[r][j];
                eign_vecs[r][i] = c*eign_vecs[r][j] + s*eign_vecs[r][i];
        }

        /* Checking for Convergence */
        off_diff = 0;
        for (r = 0; r < n; r++) {
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
        for (r = 0; r < n; r++) {
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
    for (i = 0; i < n; i++) {
        eign_vals[i] = A_tag[i][i];
    }

    /* Free A_tag <3 */
    free(A_tag[0]);
    free(A_tag);

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
    assert(values_copy != NULL && _MEM_ALLOC_ERR);
    for (i = 0; i < n; i++) {
        values_copy[i] = values[i];
    }

    /* Sort the values */   
    qsort(values, n, sizeof(int), cmpfunc);

    /* Copy the vectors after transpose into an array which allow us to swap rows (columns) in O(1) */
    vectors_T = calloc(n, sizeof(int *));
    assert(vectors_T != NULL && _MEM_ALLOC_ERR);
    for (i = 0; i < n; i++) {
        vectors_T[i] = calloc(n, sizeof(int));
        assert(vectors_T[i] != NULL && _MEM_ALLOC_ERR);
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

   printf("%f%f", vectors[0][0], values[0]); /* For compilations purposes only */
}

/* Define a comparator for the qsort function.
   Credit: 'tutorialspoint' */
int cmpfunc (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}
