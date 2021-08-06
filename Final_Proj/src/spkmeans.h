

/* Structures Declarations */
typedef struct Vector {
    int dim;
    double *values;
} Vector;

typedef struct Matrix {
    int rows;
    int cols;
    double **values;
} Matrix;

typedef struct Graph {
    double **vertices; /* Type changed from Vector. Prototype */
    double **weights;
    double *degrees;
    int n, dim;
} Graph;


/* Functions Declarations*/
double** compute_wam(Graph *g);
double compute_distance(Vector *vec1, Vector *vec2m);
double* compute_ddg(Graph *g);
double compute_degree(double **weights, int v_idx, int n);
double** compute_lnorm(Graph *g);
double* inverse_sqrt_vec(double* degs,int n);
void multi_vec_mat(double *vec, double **mat, int n, double **res);
void multi_mat_vec(double **mat, double *vec, int n, double **res);
void jacobi_alg(double **A, int n, double **eign_vecs, double *eign_vals);
void sort_by_eigen_values(double **vectors, double *values, int n);
int cmpfunc (const void * a, const void * b);
void print_matrix(double **mat, int rows, int cols); 
void read_data(char *file_path, Graph *graph);

