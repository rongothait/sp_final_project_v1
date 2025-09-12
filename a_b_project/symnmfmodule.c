 # define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"

#define CHECK_FAILURE(failure, label) \
    do { if ((failure) == 1) goto label; } while (0)

/**
 * pyListToPointList - Given a list of lists in python, it converts it to a linked list of points (and cords) in C
 * @outer_list: the python list
 * @n: the length of the python list (outer)
 * @head_point: out parameter. Will hold the first point in the linked list.
 */
static int pyListToPointList(PyObject *outer_list, int n, point **head_point){
    point *curr_point;
    cord *head_cord, *curr_cord;
    int dim, i, j;
    double num;
    PyObject *item;

    /* init head of points linked list */
    *head_point = (point*) calloc(1, sizeof(point));
    if (head_point == NULL) {return 1;}
    curr_point = *head_point;

    /* Iterating over the python list */
    for (i = 0; i < n; i++){
        PyObject *inner_list = PyList_GetItem(outer_list, i);  /* the cords for the i-th point*/
        if (!PyList_Check(inner_list)) {return 1;}

        /* init head cord for the first cord of the current point */
        head_cord = (cord*) calloc(1, sizeof(cord));
        if (head_cord == NULL) {return 1;}
        curr_cord = head_cord;

        dim = PyObject_Length(inner_list);

        for (j = 0; j < dim; j++){
            item = PyList_GetItem(inner_list, j);
            num = PyFloat_AsDouble(item);
            curr_cord->value = num;
            if (j < dim - 1){
                curr_cord->next = (cord*) calloc(1, sizeof(cord));
                if (curr_cord->next == NULL) {return 1;}
                curr_cord = curr_cord->next;
                curr_cord->next = NULL;
            }
        }
        curr_point->cords = head_cord;
        curr_point->dim = dim;

        if (i < n - 1){
            curr_point->next = (point*) calloc(1, sizeof(point));
            if (curr_point->next == NULL) {return 1;}
            curr_point = curr_point->next;
        }
    }
    return 0;  /* OK! */
}

/**
 * pyListTo2dArr - creates a 2D array from the python list of lists
 * @outer_list: the python list
 * @m: number of rows in matrix
 * @n: number of columns in matrix
 * @mat: out parameter. Will hold the C matrix result
 */
static int pyListTo2dArr(PyObject *outer_list, int m, int n, double ***mat){
    int i, j;
    PyObject *inner_list, *item;
    double num;

    /* Allocate memory for matrix */
    *mat = (double**) calloc(m, sizeof(double*));
    if (*mat == NULL) return 1;

    for (i = 0; i < m; i++){

        (*mat)[i] = (double*)calloc(n, sizeof(double));
        if ((*mat)[i] == NULL){
            /* free previously allocated memory in case of error */
            for (j = 0; j < i; j++){
                free((*mat)[j]);
            }
            free(*mat);
            return 1;
        }
    }

    /* Iterating over the python list */
    for (i = 0; i < m; i++){
        inner_list = PyList_GetItem(outer_list, i);
        if (!PyList_Check(inner_list)) { return 1;}  /* Error: inner element is not a list */

        for (j = 0; j < n; j++){
            item = PyList_GetItem(inner_list, j);
            num = PyFloat_AsDouble(item);
            if (PyErr_Occurred()) {return 1;} /* Error: element is not a float */

            (*mat)[i][j] = num;
        }
    }
    return 0;
}

static PyObject* matrix_to_pylist(double **matrix, int row_num, int col_num){
    PyObject *py_matrix, *py_row;
    int i, j;

    /* create the outer list */
    py_matrix = PyList_New(row_num);
    if (!py_matrix) return NULL;

    for (i = 0; i < row_num; i++){
        /* create the inner list */
        py_row = PyList_New(col_num);
        if (!py_row){
            Py_DECREF(py_matrix);
            return NULL; 
        }

        for (j = 0; j < col_num; j++){
            PyObject *py_value = PyFloat_FromDouble(matrix[i][j]);
            if (!py_value){
                Py_DECREF(py_row);
                Py_DECREF(py_matrix);
                return NULL;
            }
            PyList_SET_ITEM(py_row, j, py_value);
        }
        PyList_SET_ITEM(py_matrix, i, py_row);
    }

    return py_matrix;
}

static PyObject* request_standard(PyObject *args, char* goal){
    point *points_list = NULL;
    int n, failure;
    PyObject *py_points_list;
    double **ret_mat;

    /* This parses the python arguments into its C form */
    if (!PyArg_ParseTuple(args, "O", &py_points_list)) {goto error;}

    /* length of dataset */
    n = PyObject_Length(py_points_list);

    /* create linked list of points */
    failure = pyListToPointList(py_points_list, n, &points_list);
    CHECK_FAILURE(failure, error);
    
    if (strcmp(goal, "sym") == 0){
        failure = create_sim_mat(points_list, n, &ret_mat);
    }
    else if (strcmp(goal, "ddg") == 0){
        failure = create_diag_mat(points_list, n, &ret_mat);
    }
    else if (strcmp(goal, "norm") == 0){
        failure = create_normalizec_sim_mat(points_list, n, &ret_mat);
    }
    CHECK_FAILURE(failure, error);

    PyObject *py_result = matrix_to_pylist(ret_mat, n, n);
    if (!py_result) goto error;

    /* free memory */
    free_pnt_lst(points_list);
    free_matrix(ret_mat, n);

    return py_result;

error:
    free_pnt_lst(points_list);
    if (ret_mat) {free_matrix(ret_mat, n);}
    return Py_BuildValue("s", "An Error Has Occured");
}

static PyObject* sym(PyObject *self, PyObject *args){
   return request_standard(args, "sym");
}

static PyObject* norm(PyObject *self, PyObject *args){
    return request_standard(args, "norm");
}

static PyObject* ddg(PyObject *self, PyObject *args){
    return request_standard(args, "ddg");
}

static PyObject* symnmf(PyObject *self, PyObject *args){
    double **init_h_mat, **w_mat, **ret_mat;
    int n, k, failure;
    PyObject *py_w_mat, *py_init_h_mat;

    /* This parses the python arguments into its C form */
    if (!PyArg_ParseTuple(args, "OOi", &py_init_h_mat, &py_w_mat, &k)) goto error;

    n = PyObject_Length(py_w_mat);

    failure = pyListTo2dArr(py_w_mat, n, n, &w_mat);
    CHECK_FAILURE(failure, error);

    failure = pyListTo2dArr(py_init_h_mat, n, k, &init_h_mat);
    CHECK_FAILURE(failure, error);
        
    /* run the potimization algoritm of H */
    failure = optimize_h_mat(init_h_mat , w_mat , n, k, &ret_mat);
    CHECK_FAILURE(failure, error);

    /* convert and pass result to python */
    PyObject *py_result = matrix_to_pylist(ret_mat, n, k);
    if (!py_result) goto error;

    /* free memory */
    free_matrix(ret_mat, n);
    free_matrix(w_mat, n);
    free_matrix(init_h_mat, n);

    return py_result;

error:
    if (ret_mat) {free_matrix(ret_mat, n);}
    if (w_mat) {free_matrix(w_mat, n);}
    if (init_h_mat) {free_matrix(init_h_mat, n);}
    return Py_BuildValue("s", "An Error Has Occured");
}

static PyMethodDef symnmf_methods[] = {
    {
        "sym",  /* python method name */
        (PyCFunction) sym, /* C function implementing the method */
        METH_VARARGS,  /* Accepts a tuple of arguments */
        PyDoc_STR(
            "sym(dataset)\n"
            "Returns the similarity matrix of the given dataset\n"
            "parameters:\n"
            "   dataset (list of list of float): the N data points \n"
            "returns: \n"
            "   list: the resulting matrix as a list of list"
        )
    },
    {
        "ddg",  /* python method name */
        (PyCFunction) ddg, /* C function implementing the method */
        METH_VARARGS,  /* Accepts a tuple of arguments */
        PyDoc_STR(
            "ddg(dataset)\n"
            "Returns the Diagonal Degree matrix of the given dataset\n"
            "parameters:\n"
            "   dataset (list of list of float): the N data points \n"
            "returns: \n"
            "   list: the resulting matrix as a list of list"        )
    },
    {
        "norm",  /* python method name */
        (PyCFunction) norm, /* C function implementing the method */
        METH_VARARGS,  /* Accepts a tuple of arguments */
        PyDoc_STR(
            "norm(dataset)\n"
            "Returns the normalized similarity matrix of the given dataset\n"
            "parameters:\n"
            "   dataset (list of list of float): the N data points \n"
            "returns: \n"
            "   list: the resulting matrix as a list of list"        )
    },
    {
        "symnmf",  /* python method name */
        (PyCFunction) symnmf, /* C function implementing the method */
        METH_VARARGS,  /* Accepts a tuple of arguments */
        PyDoc_STR(
            "symnmf(dataset)\n"
            "Returns the optimized H matrix the given dataset\n"
            "parameters:\n"
            "   h_mat: the initial H matrix\n"
            "   w_mat: the W matrix\n"
            "   k: number of clusters\n"
            "returns:\n"
            "   list: the resulting matrix as a list of list"
        )
    },
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef symnmf_capi = {
    PyModuleDef_HEAD_INIT,
    "symnmf_mod", /* name of the module */
    NULL, /* module documentaion, may be null */
    -1, /* size of per-interpreter state of the module, -1 if the module keeps state in global variables */
    symnmf_methods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmf_mod(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmf_capi);
    if (!m){
        return NULL;
    }
    return m;
}