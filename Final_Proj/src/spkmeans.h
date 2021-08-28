/* Structures Declarations */
typedef struct Graph
{
    double **vertices;
    double **weights;
    double **lnorm;
    double *degrees;
    int N, dim;
} Graph;

typedef enum goal
{
    spk = (int)'s',
    wam = (int)'w',
    ddg = (int)'d',
    lnorm = (int)'l',
    jacobi = (int)'j'
} goal;

/* Kmeans' Structure Declaration */
struct Cluster
{
    double *vector_sum;
    int *count;
};
typedef struct Cluster Cluster;

/* SPKmeans Core Functions */
void read_data(Graph *graph, char *file_path);
void compute_by_goal(Graph *graph, goal goal);
void compute_wam(Graph *graph);
void compute_ddg(Graph *graph);
void compute_lnorm(Graph *graph);
void compute_jacobi(double **A, int N, double **eign_vecs, double *eign_vals);
void sort_by_eigen_values(double **vectors, double *values, int n);
int get_heuristic(double *eigenvalues, int N);
void form_U(double **U, double **eigenvectors, double *eigenvalues, double *eigenvalues_sorted, int N, int K);
void form_T(double **U, int N, int K);
double **init_spk_datapoints(Graph *graph, int *K);
double compute_distance_spk(double *vec1, double *vec2, int dim);
double compute_degree(double **weights, int v_idx, int n);
double compute_off_diagonal_difference(double **A, double **A_tag, int N, int i, int j);
void update_A_tag(double **A, double **A_tag, int N, int i, int j, double c, double s);
void update_A(double **A, double **A_tag, int N, int i, int j);

/* Matrix Utils */
void print_matrix(double **mat, int rows, int cols);
void print_transpose_matrix(double **mat, int rows, int cols);
void print_vector_as_matrix(double *diag, int n);
void multi_vec_mat_vec(double *vec, double **mat, int n, double **res);
void inverse_sqrt_vec(double *vector, int N, double *inv_sqrt_vec);
void init_idendity_matrix(int N, double** matrix);

/* Memory Utils */
double *calloc_vector(int size);
double **calloc_matrix(int rows, int cols);
void free_matrix(double **array);
void free_graph(Graph *graph, goal goal);

/* General Utils */
int get_sign(double d);
int cmp_func(const void *a, const void *b);
void my_assert(int status);

/* Kmeans Utils */
void kmeans(double **data_points, int N, int dim, int K, int max_iter, double **centroids);
void init_centroids(double **data_points, int K, int dim, double **centroids);
Cluster *init_clusters(int K, int dim);
int find_closest_centroid(double **centroids, double *data_point, int K, int d);
void add_datapoint_to_cluster(Cluster *clusters, int cluster_index, double *data_point, int dim);
int update_centroids(double **centroids, Cluster *clusters, int K, int dim);
double compute_distance(double *u, double *v, int dim);
