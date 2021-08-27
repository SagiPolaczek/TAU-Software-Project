#include "kmeans.c"
#include "spkmeans.h"
#include "debugger.h"
#include <math.h>

#define ERR_MSG "An Error Has Occured"

/*
    TODO:
    - Add files 'utils.c','utils.h' and move functions to them.
    - Document more.
    - Optimization!
    - Format: maximum line length: 80, maximum lines in a function: 60~
    - Check if we free memory ASAP.
    - 

    * Search 'TODO' for other tasks
    * Add 'TODO' throughout the code if there is any task :) 
*/

int MAX_ITER = 300;


/*
    Command Line Interface in the following format:
    argv = [K, goal, file_name], argc = 4
*/
int main(int argc, char *argv[]) {
    char *goal_string, *file_path;
    double **data_points, **centroids;
    int K=1, N;
    goal goal;
    Graph graph = {0};

    LOG("\n----- DEBUGGING MODE ------\n\n");
    
    if (argc != 4 && K /*compilation*/) {
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
    if (goal != spk) {
        compute_by_goal(&graph, goal);
    }
    /* goal == spk */
    else {
        data_points = init_spk_datapoints(&graph, &K);
       
        centroids = calloc_2d_array(K, K);
        init_centroids(data_points, K, K, centroids);

        kmeans(data_points, N, K, K, MAX_ITER, centroids);

        print_matrix(centroids, K, K);

        free_graph(&graph, goal);
        free_2d_array(centroids);
        free_2d_array(data_points);
    }

    return 42;
}

/*
    Read the file and store the data in the graph
*/
void read_data(Graph *graph, char *file_path) {
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
    data_points = calloc_2d_array(data_count, dim); 

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
    graph->N = data_count;    
    return;
}

/*
    Compute the Weighted Adjacency Matrix and store it in the graph.
*/
void compute_wam(Graph *graph) {
    int i, j, N, dim;
    double weight, distance;
    double *v1, *v2;
    double **vertices = graph -> vertices;
    double **weights;

    LOG("-- Computing WAM --\n");
    
    N = graph->N;
    dim = graph->dim;

    /* Allocate memory for the WAM */
    weights = calloc_2d_array(N, N);
    

    graph->weights = weights;
    
    for (i = 0; i < N; i++) {
        v1 = vertices[i];

        for (j = 0; j < i; j++) {
            v2 = vertices[j];
            distance = compute_distance_spk(v1, v2, dim);
            weight = exp((-1) * (distance / 2));
            weights[i][j] = weight;
            weights[j][i] = weight;
        }
        weights[i][i] = 0;
    }
    LOG("-- COMPLETE COMPUTE WAM --\n");
    return;
}


/*
    Compute the distance between the two vector.
    The vectors have 'dim' as dimention.
 */
double compute_distance_spk(double *vec1, double *vec2, int dim) {
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
void compute_ddg(Graph *graph) {
    int i, N;
    double deg;
    double **weights, *ddg;

    LOG("-- Computing DDG --\n");
    weights = graph -> weights;
    N = graph -> N;

    /* Allocate memory for the DDG's "matrix" - represented by vector */
    ddg = calloc_1d_array(N);
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
double compute_degree(double **weights, int v_idx, int n) {
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
void compute_lnorm(Graph *graph) {
    double *invsqrt_d;
    int i, j;
    double **weights = graph -> weights;
    double *degs = graph -> degrees;
    int N = graph -> N;
    double **lnorm;

    LOG("-- Computing LNORM --\n");

    /* Allocate memory for the Lnorm's matrix */
    lnorm = calloc_2d_array(N, N);
    graph->lnorm = lnorm;

    invsqrt_d = calloc_1d_array(N);

    inverse_sqrt_vec(degs, N, invsqrt_d);

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

void inverse_sqrt_vec(double *vector, int N, double *inv_sqrt_vec) {
    int i;
    
    for (i = 0; i < N; i++) {
        inv_sqrt_vec[i] = 1 / (sqrt(vector[i]));
    }
    return;
}

/*
    Compute multiplication of: ' VECTOR * MATRIX '
    when the 'VECTOR' is a diagonal matrix A such that A[i][i]=V[i]
*/
void multi_vec_mat(double *vec, double **mat, int n, double **res) {
    int i, j;

    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            res[i][j] = vec[i] * mat[i][j];
        }
    }
}

/*
    Compute multiplication of: ' MATRIX * VECTOR '
    when the 'VECTOR' is a diagonal matrix A such that A[i][i]=V[i]
*/
void multi_mat_vec(double **mat, double *vec, int n, double **res) {
    int i, j;

    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            res[j][i] = vec[i] * mat[j][i];
        }
    } 
}

/*
    TODO:
    Elaborate & split to several functions (sagi).
*/
void compute_jacobi(double **A, int N, double **eigen_vecs, double *eigen_vals) {
    int is_not_diag = 1;
    int i, j, r;
    double max_entry, curr_entry;
    int max_row, max_col;
    int sign_theta;
    double theta, c, s, t, s_sq, c_sq;
    double **A_tag, **eigen_vecs_dcopy;
    double off_diff;
    int iter_count;
    const int MAX_ITER = 100;
    const double EPS = 0.001;

    LOG("-- Computing JACOBI --\n");

    /* A_tag start as deep copy of A */
    A_tag = calloc_2d_array(N, N);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A_tag[i][j] = A[i][j];
        }
    }
    
    /* Initalize the Idendity Matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            eigen_vecs[i][j] = (i == j);
        }
    }

    /* Eigenvaluse deep copy */
    eigen_vecs_dcopy = calloc_2d_array(N, N);
    
    
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
        theta = (A[j][j] - A[i][i]) / (2*A[i][j]);
        sign_theta = get_sign(theta);
        t = sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / sqrt(pow(t, 2) + 1);
        s = t * c;
        
        /* Relation between A and A_tag */
        i = max_row;
        j = max_col;
        for (r = 0; r < N; r++) {
            if ((r != i) && (r != j)) {
                A_tag[r][i] = c*A[r][i] - s*A[r][j];
                A_tag[r][j] = c*A[r][j] + s*A[r][i];
                A_tag[i][r] = A_tag[r][i];
                A_tag[j][r] = A_tag[r][j];
            }
        }
        c_sq = pow(c, 2);
        s_sq = pow(s, 2);
        A_tag[i][i] = c_sq*A[i][i] + s_sq*A[j][j] - 2*s*c*A[i][j];
        A_tag[j][j] = s_sq*A[i][i] + c_sq*A[j][j] + 2*s*c*A[i][j];
        A_tag[i][j] = 0; /* (c_sq - s_sq)*A[i][j] + s*c*(A[i][i] - A[j][j]) => =0 */
        A_tag[j][i] = 0;

        /* Update the deep copy of the eigenvectors*/
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
            eigen_vecs_dcopy[i][j] = eigen_vecs[i][j];
            }
        }

        i = max_row;
        j = max_col;
        /* Update Eigenvectors' matrix */
        for (r = 0; r < N; r++) {
                eigen_vecs[r][i] = c*eigen_vecs_dcopy[r][i] - s*eigen_vecs_dcopy[r][j];
                eigen_vecs[r][j] = c*eigen_vecs_dcopy[r][j] + s*eigen_vecs_dcopy[r][i];
        }

        /* Checking for Convergence */
        off_diff = 0;
        for (r = 0; r < N; r++) {
            off_diff += pow(A[r][i], 2);
            off_diff += pow(A[r][j], 2);
            off_diff -= pow(A_tag[r][i], 2);
            off_diff -= pow(A_tag[r][j], 2);
        }
        off_diff -= pow(A[i][i], 2);
        off_diff -= pow(A[j][j], 2);
        off_diff += pow(A_tag[i][i], 2);
        off_diff += pow(A_tag[j][j], 2);

        if (off_diff <= EPS || iter_count > MAX_ITER) {
            is_not_diag = 0;
            /* 'Deep' update A = A' */
            for (r = 0; r < N; r++) {
                A[r][i] = A_tag[r][i];
                A[r][j] = A_tag[r][j];
                A[i][r] = A_tag[r][i];
                A[j][r] = A_tag[r][j];
            }
            break;
        }
        
        /* 'Deep' update A = A' */
        for (r = 0; r < N; r++) {
            A[r][i] = A_tag[r][i];
            A[r][j] = A_tag[r][j];
            A[i][r] = A_tag[r][i];
            A[j][r] = A_tag[r][j];
        }

        iter_count += 1;
    }

    /* Update Eigenvalues */
    for (i = 0; i < N; i++) {
        eigen_vals[i] = A[i][i];
    }

    /* Free <3 */
    free_2d_array(A_tag);
    free_2d_array(eigen_vecs_dcopy);
    return;
}

/*
    Define a comperator for the qsort function.
    SOURCE: 'https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm'
*/
int cmp_func (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

/*
    Print the matrix with the given dimention in the right format.
*/
void print_matrix(double **mat, int rows, int cols) {
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
        printf("\n");
    }
}

/*
    Print the matrix with the given dimention in the right format,
    when PRINTED[i][j] = MATRIX[j][i].

    NOTE:
    We define a special function for that because we didn't want to spend time on rotating the matrix.
    Its a tradeoff of: Time Vs Additional code and elegant. We chose time :) .
*/
void print_transpose_matrix(double **mat, int rows, int cols) {
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
void print_vector_as_matrix(double *diag, int n) {
    double val;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                val = diag[i];
            } else {
                val = 0;
            }
            printf("%.4f", val);
            if (j < n - 1) {
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
double **init_spk_datapoints(Graph *graph, int *K) {
    double **eigenvectors, *eigenvalues, *eigenvalues_sorted;
    int N, i;
    double **U, **T;

    /* FLOW:
    1) call compute wam -> ddg -> lnorm -> jacobi -> (?) U -> T
    2) return T
    */
    compute_wam(graph);
    compute_ddg(graph);
    compute_lnorm(graph);

    /* Allocate memory for the eigenvectors & eigenvalues */
    N = graph->N;
    eigenvectors = calloc_2d_array(N, N);
    eigenvalues = calloc_1d_array(N);

    compute_jacobi(graph->lnorm, N, eigenvectors, eigenvalues);

    /* Deep copy the eigenvalues */
    eigenvalues_sorted = calloc_1d_array(N);

    for (i = 0; i < N; i++) {
        eigenvalues_sorted[i] = eigenvalues[i];
    }

    /* Sort eigenvalues_sorted */
    qsort(eigenvalues_sorted, N, sizeof(int), cmp_func);

    /*  If K doesn't provided (k==0) determine k */
    if (K == 0) {
        *K = get_heuristic(eigenvalues, N);
    }

    /*  Obtain the first (ordered!) k eigenvectors u1, ..., uk
        Form the matrix U(nxk) which u_i is the i'th column */

    /* Allocate memory for the U matrix. nxk. */
    U = calloc_2d_array(N, (*K));

    form_U(U, eigenvectors, eigenvalues, eigenvalues_sorted, N, *K);     

    /*  Form T from U by renormalizing each row to have the unit length */
    T = U;
    form_T(T, *K, N);

    free(eigenvalues_sorted);
    free(eigenvalues); 
    free_2d_array(eigenvectors);
    
    return T;
}


/*
    TODO:
    elaborate.
*/
int get_heuristic(double *eigenvalues, int N) {
    double max_delta;
    int max_idx = ((N/2)+1 < N-2)? (N/2)+1 : N-2; /* i is bounded as shown in the pdf */
    int i, K;
    double val1, val2, curr_delta;

    for (i = 0; i < max_idx; i++) {
        val1 = eigenvalues[i];
        val2 = eigenvalues[i+1];
        curr_delta = val1 - val2;
        if (curr_delta > max_delta) {
            max_delta = curr_delta;
            K = i;
        }
    }
    return K;
}

/*
    TODO:
    elaborate.
*/
void form_U(double **U, double **eigenvectors, double *eigenvalues, double *eigenvalues_sorted, int N, int K) {
    double curr_val, curr_sorted_val;
    const int NULL_VAL = -42; /* define null value when known that eigenvalues are positive */
    int i, j, r;

    for (i = 0; i < K; i++) {
        curr_sorted_val = eigenvalues_sorted[i];
        for (j = 0; j < N; j++) {
            curr_val = eigenvalues[j];
            if (curr_val == curr_sorted_val) {
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
    TODO:
    elaborate.
*/
void form_T(double **U, int N, int K) {
    int i, j;
    double vec_len;
    double *vec, *zero_vec;

    /* init the zero_vec */
    zero_vec = calloc_1d_array(K);

    for (i = 0; i < N; i++) {
        vec = U[i];
        /* compute the row vector's length */
        vec_len = compute_distance_spk(vec, zero_vec, K);

        for (j = 0; j < K; j++) {
            vec[j] = ( vec[j] / vec_len );
        }
    }

    free(zero_vec);
}

/*
    Allocating memory for a vector.
    TODO:
    change name to 'calloc_vector'.
*/
double *calloc_1d_array(int size) {
    double *array;
    array = calloc(size, sizeof(double));
    my_assert(array != NULL);

    return array;
}

/*
    Allocating memory for a matrix.
    TODO:
    change name to 'calloc_matrix'.
*/
double **calloc_2d_array(int rows, int cols) {
    double *ptr, **array;
    int i;

    ptr = calloc((rows * cols), sizeof(double));
    array = calloc(rows, sizeof(double*));
    my_assert(array != NULL || ptr != NULL);
    for (i = 0; i < rows; i++) {
        array[i] = ptr + i*cols;
    }

    return array;
}

/*
    Free a memory of a matrix.
    TODO:
    change name to 'free_matrix'.
*/
void free_2d_array(double **array) {
    free(array[0]);
    free(array);
}

/*
    Free the memory of the graph.
*/
void free_graph(Graph *graph, goal goal) {
    LOG("-- FREE GRAPH --\n");

    free_2d_array(graph->vertices);
    if (goal == jacobi) {return;}

    free_2d_array(graph->weights);
    printf("TEST TEST %d\n", (goal == wam));
    if (goal == wam) {return;}

    printf("TEST TEST\n");
    free(graph->degrees);

    if (goal == ddg) {return;}

    free_2d_array(graph->lnorm);
    LOG("HERE\n");

    return;
}

/*
    Return the sign of a double.
*/
int get_sign(double d) {
    if (d >= 0)
        return 1;
    else
        return -1;
}

/*
    Custom assert fucntion to satisfy the assignment's requirements.
*/
void my_assert(int status) {
    if (status == 0) {
        printf(ERR_MSG);
        exit(1);
    }
}

void compute_by_goal(Graph *graph, goal goal) {
    double *eigenvalues, **eigenvectors, **A;
    int N = graph->N;
    if (goal == jacobi) { 

        /* Allocate memory for the eigenvectors & eigenvalues */
        eigenvectors = calloc_2d_array(N, N);
        eigenvalues = calloc_1d_array(N);

        A = graph->vertices;
        compute_jacobi(A, N, eigenvectors, eigenvalues);

        /* Print eigenvectors & eigenvalues */
        LOG("-- Eigenvalues:\n");
        print_matrix(&eigenvalues, 1, N);

        LOG("-- Eigenvectors:\n");
        print_transpose_matrix(eigenvectors, N, N);

        free_2d_array(eigenvectors);
        free(eigenvalues);

        free_graph(graph, jacobi);

        LOG("-- COMPLETE GOAL JACOBI --");
        return;
    }

    compute_wam(graph);

    /* If the goal is only to compute the WAM, we finished */
    if (goal == wam) {
        print_matrix(graph->weights, N, N);

        free_graph(graph, wam);

        printf("C - %d\n", 10);
        LOG("-- COMPLETE GOAL WAM --");
        return;
    }

    compute_ddg(graph);
    
    if (goal == ddg) {
        print_vector_as_matrix(graph->degrees, N);

        free_graph(graph, ddg);
        LOG("-- COMPLETE GOAL DDG --");
        return;
    }

    compute_lnorm(graph);

    if (goal == lnorm) {
        print_matrix(graph->lnorm, N, N);

        free_graph(graph, lnorm);
        LOG("-- COMPLETE GOAL LNORM --");
        return;
    }
}
