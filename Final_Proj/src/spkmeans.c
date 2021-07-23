#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>


/* Structures Declarations */
struct Vector {
    int dim;
    double *values;
} typedef Vector;

struct Matrix {
    int rows;
    int cols;
    double **values;
} typedef Matrix;

struct Graph {
    Vector **vertices;
    double **weights;
    double *degrees;
    int n, dim;
} typedef Graph;





/* Functions Declarations */
double** compute_spk();
double** compute_wam(Graph *g);
double compute_distance(Vector *vec1, Vector *vec2m);
double* compute_ddg(Graph *g);
double compute_degree(double **weights, int v_idx, int n);
double** compute_lnorm(Graph *g);
double** compute_jacobi();
static PyObject* fit(PyObject* self, PyObject* args);


/* Python */


/* C */
double **compute_spk() {
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
    int i, j;
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

void muilti_mat_vec(double **mat, double *vec, int n, double **res) {
    int i, j;

    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            res[j][i] = vec[i] * mat[j][i];
        }
    } 
}