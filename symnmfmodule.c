 # define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"

#define CHECK_FAILURE(failure, label) \
    do { if ((failure) == 1) goto label; } while (0)

static const char *ERR_MSG_GENERAL = "An Error Has Occurred";

/**
 * calculates the length of every inner list in a 2D py list
 * @outer_list: A 2D python lists of lists
 * @n: number of lists in the list
 * @dim: out parameter. will hold the length of the inner lists
 */
static int calc_cols(PyObject *outer_list, int n, int *dim){
    int expected_cols, i, curr_col;
    PyObject *inner_list;
    expected_cols = -1;

    for (i = 0; i < n; i++){
        inner_list = PyList_GetItem(outer_list, i);
        if (!PyList_Check(inner_list)) { goto error; }
        curr_col = PyObject_Length(inner_list);
        if (expected_cols == -1) expected_cols = curr_col;
        else if (expected_cols != curr_col) goto error; /* dimension are not the same for all of the points */
    }
    *dim = expected_cols;

error:
    return -1;
}

/**
 * convert python list of lists to 2d array
 * @outer_list: python list of lists
 * @n: number of lists
 * @cols_out: out_parameter. will hold the dimension of every list.
 * @dataset_mat: out parameter. will hold the matrix with the dataset points.
 */
static int py_list_to_matrix(PyObject *outer_list, int n, int *cols_out, double ***dataset_mat){
    int i,j, failure;
    double num;
    PyObject *item, *inner_list;
    
    failure = calc_cols(outer_list, n, cols_out);
    CHECK_FAILURE(failure, error);

    failure = allocate_double_matrix(n, *cols_out, dataset_mat);
    CHECK_FAILURE(failure, error);

    for (i = 0; i < n; i++){
        inner_list = PyList_GetItem(outer_list, i);
        if (!PyList_Check(inner_list)) goto error; 
        for (j = 0; j < *cols_out; j++){
            item = PyList_GetItem(inner_list, j);
            num = PyFloat_AsDouble(item);
            (*dataset_mat)[i][j] = num;
        }
    }
    return 0;

error:
    return 1;
}

/**
 * matrix_to_pylist - converts a C 2D array (matrix) to a python list of lists
 * @matrix: the input C 2d array (double**)
 * @row_num: Number of rows in the matrix
 * @col_num: number of columns in the matrix
 */
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

/**
 * request_standard - Handles standard matrix goals (sym, norm, ddg) for the Python extension
 * @args: Python tuple containing the dataset as a list of lists.
 * @goal: string specifiying the matrix operation ("sym", "ddg", "norm")
 */
static PyObject* request_standard(PyObject *args, char* goal){
    double **dataset_mat = NULL;
    int n = -1, failure, dimension;
    PyObject *py_dataset;
    double **ret_mat = NULL;

    if (!PyArg_ParseTuple(args, "O", &py_dataset)) { goto error; } /* This parses the python arguments into its C form */

    n = PyObject_Length(py_dataset); /* length of dataset */
    failure = py_list_to_matrix(py_dataset, n, &dimension, &dataset_mat);  /* create matrix of dataset points */
    CHECK_FAILURE(failure, error);
    
    if (strcmp(goal, "sym") == 0){
        failure = create_sim_mat(dataset_mat, n, dimension, &ret_mat);
    }
    else if (strcmp(goal, "ddg") == 0){
        failure = create_diag_mat(dataset_mat, n, dimension, &ret_mat);
    }
    else if (strcmp(goal, "norm") == 0){
        failure = create_normalized_sim_mat(dataset_mat, n, dimension, &ret_mat);
    }
    CHECK_FAILURE(failure, error);

    PyObject *py_result = matrix_to_pylist(ret_mat, n, n);
    if (!py_result) goto error;

    free_matrix(&dataset_mat, n); /* free memory */
    free_matrix(&ret_mat, n);
    return py_result;

error:
    if (dataset_mat) free_matrix(&dataset_mat, n);
    if (ret_mat) {free_matrix(&ret_mat, n);}
    PyErr_SetString(PyExc_ValueError, ERR_MSG_GENERAL);
    return NULL;
}

/**
 * sym - Python exposed function to compute for the sym call
 */
static PyObject* sym(PyObject *self, PyObject *args){
   return request_standard(args, "sym");
}

/**
 * norm - Python exposed function to compute for the norm call
 */
static PyObject* norm(PyObject *self, PyObject *args){
    return request_standard(args, "norm");
}

/**
 * ddg - Python exposed function to compute for the ddg call
 */
static PyObject* ddg(PyObject *self, PyObject *args){
    return request_standard(args, "ddg");
}

/**
 * symnmf - Python exposed function to compute for the symnmf call
 */
static PyObject* symnmf(PyObject *self, PyObject *args){
    double **init_h_mat = NULL, **w_mat = NULL, **ret_mat = NULL;
    int n = -1, k, failure, w_cols;
    PyObject *py_w_mat, *py_init_h_mat;

    /* This parses the python arguments into its C form */
    if (!PyArg_ParseTuple(args, "OO", &py_init_h_mat, &py_w_mat)) goto error;

    n = PyObject_Length(py_w_mat);

    failure = py_list_to_matrix(py_w_mat, n, &w_cols, &w_mat);
    CHECK_FAILURE(failure, error);

    failure = py_list_to_matrix(py_init_h_mat, n, &k, &init_h_mat);
    CHECK_FAILURE(failure, error);
        
    /* run the potimization algoritm of H */
    failure = optimize_h_mat(init_h_mat , w_mat , n, k, &ret_mat);
    CHECK_FAILURE(failure, error);

    /* convert and pass result to python */
    PyObject *py_result = matrix_to_pylist(ret_mat, n, k);
    if (!py_result) goto error;

    /* free memory */
    free_matrix(&ret_mat, n);
    free_matrix(&w_mat, n);
    free_matrix(&init_h_mat, n);

    return py_result;

error:
    if (ret_mat) {free_matrix(&ret_mat, n);}
    if (w_mat) {free_matrix(&w_mat, n);}
    if (init_h_mat) {free_matrix(&init_h_mat, n);}
    PyErr_SetString(PyExc_ValueError, ERR_MSG_GENERAL);
    return NULL;
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