#define PY_SSIZE_T_CLEAN
#include <Python.h>

int N, k, d, max_iter;

static double ** getInput(PyObject * vectors);
static double dist(double * vec1, double * vec2);
static double ** kmeans(PyObject * observations, PyObject * firsts);
static int belongs (double ** centroids, double * vec);
static double * addVectors(double * vec1, double * vec2);
static double * centroidCalc(double * vec, int member);
static PyObject * kmeanspp(PyObject *self, PyObject *args);
static void print_centroids(double ** centroids);
static void print_indexes(PyObject * list, Py_ssize_t list_size);



static PyObject * kmeanspp(PyObject *self, PyObject *args)
{
    PyObject *vectors, *centroids, *vector1, *vector2, *indexes;
    double ** output;
    Py_ssize_t i;
    if(!PyArg_ParseTuple(args, "OOOiiii:Problem with accepting arguments from Python", &vectors, &centroids, &indexes, &max_iter, &k, &N, &d)) {
        PyErr_SetString(PyExc_Exception, "There was a problem with the input!");
        return NULL;
    }
    /* Is it a list? */
    if (!PyList_Check(vectors) || !PyList_Check(centroids)){
        PyErr_SetString(PyExc_Exception, "There was a problem with the input!");
        return NULL;
    }
    /* Get the size of it and build the output list */
    if(!PyList_Check(indexes)){
        PyErr_SetString(PyExc_Exception, "There was a problem with the indexes");
        return NULL;
    }
    for (i = 0; i < PyList_Size(vectors); i++) {
        vector1 = PyList_GetItem(vectors, i);
        if (i < PyList_Size(centroids)) {
            vector2 = PyList_GetItem(vectors, i);
            if (!PyList_Check(vector2)) {
                PyErr_SetString(PyExc_Exception, "There was a problem with the input - it's not a list of lists!");
                return NULL;
            }
        }
        if (!PyList_Check(vector1)) {  /* We only print lists */
                PyErr_SetString(PyExc_Exception, "There was a problem with the input - it's not a list of lists!");
                return NULL;
            }
    }
    output = kmeans(vectors, centroids);
    print_indexes(indexes, PyList_Size(indexes));
    print_centroids(output);
    for (i = 0; i < PyList_Size(centroids); i++) {
        free(output[i]);
    }
    free(output);
    Py_RETURN_NONE;
}

static void print_indexes(PyObject * list, Py_ssize_t list_size){
    Py_ssize_t i;
    PyObject *item;
    long * my_c_list = malloc(sizeof my_c_list * list_size);
    assert(my_c_list != NULL && "Problem in print_list_of_ints()");
    for (i = 0; i < list_size; i++) {
        item = PyList_GetItem(list, i);
        assert(PyLong_Check(item)!=NULL);
        my_c_list[i] = PyLong_AsLong(item);
        if (my_c_list[i]  == -1 && PyErr_Occurred()){
            puts("Integer too big to fit in a C long, bail out");
            free(my_c_list);
            return;
         }
        if (i != list_size - 1) {
            printf("%ld,", my_c_list[i]);
        }
        else {
            printf("%ld", my_c_list[i]);
        }
    }
    printf("\n");
    free(my_c_list);
    return;
}

static void print_centroids(double ** centroids) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            printf("%.2f", centroids[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }

}
static double dist(double * vec1, double * vec2){
    double counter;
    int i;
    counter = 0;
    for(i = 0; i < d; i++){
        counter += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return counter;
}

static int belongs(double ** centroids, double * vec){
    int minInd, i;
    double minDist, currDist;
    minInd = 0;
    minDist = dist(centroids[0],vec);
    for (i = 1; i < k; i++){
        currDist = dist(centroids[i], vec);
        if (currDist < minDist){
            minInd= i;
            minDist = currDist;
        }

    }
    return minInd;
}

static double * addVectors(double * vec1, double * vec2){
    int i;
    for(i = 0; i < d; i++){
        vec1[i] += vec2[i];
    }
    return vec1;

}


static double * centroidCalc(double * vec, int member){
    int i;
    for(i = 0; i < d; i++){
        vec[i] = (vec[i] / member);
    }
    return vec;
}


static double ** getInput(PyObject * vectors) {
    double ** newVecs = (double **) calloc (PyList_Size(vectors), sizeof(double *));
    assert(newVecs != NULL);
    PyObject * item;
    PyObject * inside;
    Py_ssize_t i, j;
    for (i = 0; i < PyList_Size(vectors); i++) {
        item = PyList_GetItem(vectors, i);
        newVecs[i] = (double *) calloc(PyList_Size(item), sizeof(double));
        for (j = 0; j < PyList_Size(item); j++) {
            assert(newVecs[j] != NULL);
            inside = PyList_GetItem(item, j);
            if (!PyFloat_Check(inside)) {
                // free here
                printf("Problem!");
                return NULL;
            }
            newVecs[i][j] = PyFloat_AsDouble(inside);
        }
    }
    return newVecs;
}

static double ** kmeans(PyObject * observations, PyObject * firsts) {
    int ind, changed, g, n, m, y, i, j, f, p;
    double ** vectors = getInput(observations);
    double ** initials = getInput(firsts);
    int * members;
    double ** sumVectors;
    ind = 0;
    changed = 0;
    members = (int *) calloc(k, sizeof(int));
    assert(members != NULL);
    sumVectors = (double **) malloc(k * sizeof(double *));
    assert(sumVectors != NULL);
    for (g = 0; g < k; g++) {
        sumVectors[g] = (double *) calloc(d, sizeof(double));
        assert(sumVectors[g] != NULL);
    }
    while(ind < max_iter && changed == 0) {
        for (n = 0; n < N; n++) {
            int belong = belongs(initials, vectors[n]);
            members[belong] += 1;
            sumVectors[belong] = addVectors(sumVectors[belong], vectors[n]);
        }
        changed=1;
        for (m = 0; m < k; m++) {
            double *currCent = centroidCalc(sumVectors[m], members[m]);
            for (f = 0; f < d; f++) {
                if (currCent[f] != initials[m][f]) {
                    changed = 0;
                    break;
                }
            }
            for (y = 0; y < d; y++) {
                initials[m][y] = currCent[y];
            }
        }
        for (j = 0; j < k; j++) {
            members[j] = 0;
            for(i = 0; i < d; i++){
                sumVectors[j][i] = 0;
            }
        }
        ind += 1 ;
    }
    free(members);
    for (p = 0; p < N; p++) {
        free(vectors[p]);
        if (p < k) {
            free(sumVectors[p]);
        }
    }
    free(sumVectors);
    free(vectors);
    return initials;
}

static PyMethodDef myMethods[] = {
        {"kmeanspp", (PyCFunction) kmeanspp, METH_VARARGS, PyDoc_STR("Kmeanspp function")},
        {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        myMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    return PyModule_Create(&moduledef);
}