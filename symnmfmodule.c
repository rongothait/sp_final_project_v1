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

static PyObject* sym(PyObject *self, PyObject *args){
    point *points_list = NULL;
    int n, failure;
    PyObject *py_points_list;
    double **sym_mat;

    /* This parses the python arguments into its C form */
    if (!PyArg_ParseTuple(args, "O", &py_points_list)) {goto error;}

    /* length of dataset */
    n = PyObject_Length(py_points_list);

    /* create linked list of points */
    failure = pyListToPointList(py_points_list, n, &points_list);
    CHECK_FAILURE(failure, error);

    char *result_str = NULL;
    failure = create_sim_mat(points_list, n, &sym_mat);
    CHECK_FAILURE(failure, error);

    failure = sqr_mat_to_str(sym_mat, n, &result_str);
    CHECK_FAILURE(failure, error);
    PyObject *py_result = Py_BuildValue("s", result_str);

    /* free memory */
    free(result_str);
    free_pnt_lst(points_list);

    return py_result;

error:
    if (result_str) {free(result_str);}
    free_pnt_lst(points_list);
    return Py_BuildValue("s", "An Error Has Occured");
}

static PyObject* ddg(PyObject *self, PyObject *args){
    point *points_list = NULL;
    int n, failure;
    PyObject *py_points_list;
    double **ddg_mat;

    /* This parses the python arguments into its C form */
    if (!PyArg_ParseTuple(args, "O", &py_points_list)) {goto error;}

    /* length of dataset */
    n = PyObject_Length(py_points_list);

    /* create linked list of points */
    failure = pyListToPointList(py_points_list, n, &points_list);
    CHECK_FAILURE(failure, error);

    char *result_str = NULL;
    failure = create_diag_mat(points_list, n, &ddg_mat);
    CHECK_FAILURE(failure, error);

    failure = sqr_mat_to_str(ddg_mat, n, &result_str);
    CHECK_FAILURE(failure, error);
    PyObject *py_result = Py_BuildValue("s", result_str);

    /* free memory */
    free(result_str);
    free_pnt_lst(points_list);

    return py_result;

error:
    if (result_str) {free(result_str);}
    free_pnt_lst(points_list);
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
            "   str: the resulting matrix as a formatted string"
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
            "   str: the resulting matrix as a formatted string"
        )
    },
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef symnmf = {
    PyModuleDef_HEAD_INIT,
    "symnmf", /* name of the module */
    NULL, /* module documentaion, may be null */
    -1, /* size of per-interpreter state of the module, -1 if the module keeps state in global variables */
    symnmf_methods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmf);
    if (!m){
        return NULL;
    }
    return m;
}