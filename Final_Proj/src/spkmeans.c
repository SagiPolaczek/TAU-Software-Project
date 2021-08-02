#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "spkmeans.h"


int main(int argc, char *argv[]) {


    printf("%d%s", argc, argv[0]);
    return 0;
}

double **compute_wam(Graph *g) {
    int i, j, n;
    double weight, distance;
    Vector *v1, *v2;
    Vector **vertices = g -> vertices;
    double **weights = g -> weights;

    n = g -> n;
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
    }
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
    assert(p != NULL);
    res = calloc(n, sizeof(double*));
    assert(res != NULL);

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
    assert(res != NULL);
    
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
    assert(p != NULL);
    A_tag = calloc(n, sizeof(double*));
    assert(A_tag != NULL);

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
void sort_by_eigen_values(double **vectors, double *values) {
    /* Idea:
        qsort the eigenvalues. (fastest & easiest way).
        Then order the eigenvectors to correspond the eigenvalues order.
        The last part will take O(n^2) but maybe should be pretty fast.

        NOTE:
        We should change the matrix reprasentation. 
        If so, we can achive O(1) time for swapping between two vectors instead of O(n).
    */
}

