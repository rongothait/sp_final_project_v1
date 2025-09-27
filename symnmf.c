/**
 * Symnmf clustering algorithm
 * All functions return int: 0 = success, 1 = failure.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CHECK_FAILURE(failure, label) \
    do { if ((failure) == 1) goto label; } while (0)

static const char *ERR_MSG_GENERAL = "An Error Has Occurred";
int failure;
double EPSILON = 0.0001;
int MAX_ITER = 300;

/**
 * allocate_double_matrix - allocate rows x cols matrix (double**)
 * @rows: number of rows
 * @cols: number of columns
 * @ret_mat: out parameter; set to the newly allocated matrix on success
 */
int allocate_double_matrix(int rows, int cols, double ***ret_mat){
    int i,k;
    *ret_mat = NULL;

    *ret_mat = (double**)calloc(rows, sizeof(double*));
    if (*ret_mat == NULL) { return 1; }

    for (i = 0; i < rows; i++){
        (*ret_mat)[i] = (double*)calloc(cols, sizeof(double));
        if ((*ret_mat)[i] == NULL){
            /* free previousley allocated memory */
            for (k = 0; k < i; k++) { free((*ret_mat)[k]); }
            free(*ret_mat);
            *ret_mat = NULL;
            return 1;
        }
    }
    return 0;
}


/**
 * free_matrix - Frees the memory allocated for a 2d matrix mat with m rows
 * @mat: the matrix (a 2d array)
 * @m: number of rows
 */
int free_matrix(double ***pmat, int m){
    int i;
    double **mat;

    if (!pmat || !*pmat) return 0;  /* nothing to free */

    mat = *pmat;

    /* Free each row */
    for (i = 0; i < m; i++){
        free(mat[i]);
        mat[i] = NULL;
    }

    free(mat);
    *pmat = NULL; /* prevent double-free via same variable */
    return 0;
}

/**
 * mat_to_str - returns a string describing the matrix
 * @mat: n*n matrix to be descibed
 * @m: rows dimensions of matrix
 * @n: columns dimension of matrix
 * @ret_str: outparameter. Will hold the describing string
 *  Return: 0 on success, 1 on failure
 */
int mat_to_str(double** mat, int m, int n, char **ret_str){
    int i, j, len = 0;
    int bufsize = m * n * 32 + n; /* Estimate: 32 chars per number + newlines */
    *ret_str = malloc(bufsize);
    if (*ret_str == NULL) {return 1;}

    (*ret_str)[0] = '\0'; /* Start with empty string */

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            /* Append number, add comma if not last in row */
            len += sprintf(*ret_str + len, "%.4f", mat[i][j]);
            if (j < n - 1)
                len += sprintf(*ret_str + len, ",");
        }
        if (i < m - 1)
            len += sprintf(*ret_str + len, "\n");
    }
    return 0;
}

/**
 * read file and count how many rows and how many columns
 * @path: path to dataset file
 * @rows_out: out parameter. will store how many rows in dataset
 * @cols_out: out parameter. will store how many columns in dataset
 */
int calc_file_dimension(char *path, int *rows_out, int *cols_out){
    FILE *file;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int i, lines_count = 0, expected_cols = -1, curr_cols_count;

    file = fopen(path, "r");
    if (file == NULL) goto error;

    while ((read = getline(&line, &len, file)) != -1){
        if (read == 1 && line[0] == '\n') continue; /* last line is expected to be blank */

        lines_count++; /* increment number of lines in file */
        curr_cols_count = 1; /* count columns in line */
        for (i = 0; i < read; i++){
            if (line[i] == ',') curr_cols_count++;
        }
        if (expected_cols == -1) {  /* if first line set the expected cols */
            expected_cols = curr_cols_count;
        } else if (curr_cols_count != expected_cols) goto error; /* inconsistent number of columns in file - error */
    }

    *rows_out = lines_count; 
    *cols_out = expected_cols;
    
    free(line);
    fclose(file);
    return 0;

error:
    if (line) free(line);
    if (file) fclose(file);
    return 1;
}

/**
 * convert txt file into a 2d double array
 * @path: the path to the file containing the dataset
 * @rows_out: out parameter. will hold the number of rows in the dataset
 * @cols_out: out parameter. will hold the number of columns in the dataset
 * @dataset_mat: out parameter. will hild the 2D array with the data from the file
 */
int input_txt_to_matrix(char *path, int *rows_out, int *cols_out, double ***dataset_mat){
    int failure, i, j;
    FILE *file;
    char c;
    double n;
    *rows_out = -1;
    *dataset_mat = NULL;
    file = NULL;

    /* calucalte the dimension and allocate a matrix of that size*/
    failure = calc_file_dimension(path, rows_out, cols_out);
    CHECK_FAILURE(failure, error);
    failure = allocate_double_matrix(*rows_out, *cols_out, dataset_mat);
    CHECK_FAILURE(failure, error);

    file = fopen(path, "r");
    if (file == NULL) goto error;

    /* populate matrix with data */
    for (i = 0; i < *rows_out; i++){
        for (j = 0; j < *cols_out; j++){
            if (fscanf(file, "%lf%c", &n, &c) != 2) goto error;
            (*dataset_mat)[i][j] = n;
        }
    }
    
    fclose(file);
    return 0;
error:
    if (file) fclose(file);
    if (*rows_out != -1) free_matrix(dataset_mat, *rows_out);
    return 1;
}

/**
 * euc_distance - Calculates the eucledian distance between p1 and p2
 * @p1: the first point
 * @p2: the second point
 * @dimension: dimension of p1 and p2 (needs to be same dimension)
 * @dist: output parameter. Will hold the eucledian distance of the 2 points.
 *  Return: 0 on success, 1 on failure
 */
int euc_distance(double *p1, double *p2, int dimension, double *dist){
    int i;
    *dist = 0;

    for (i = 0; i < dimension; i++){
        (*dist) += pow((p1[i] - p2[i]), 2);
    }

    *dist = sqrt(*dist);
    return 0;
}

/**
 * calc_aij_sim - Calculates the value of A(ij) in the similiarity matrix as defined in 1.1
 * @pre: xi.dim == x2.dim
 * @xi: point i in dataset X
 * @xj: point j in dataset X
 * @dimension: dimension of the points
 * @result: out parameter. Will hold the result of A(ij)
 *  Return: 0 on success, 1 on failure
 */
int calc_aij_sim(double *xi, double *xj, int dimension, double *result){
    double dist, power;
    *result = 0;
    if(euc_distance(xi, xj, dimension, &dist) != 0) { return 1; } /* error! */
    power = ((pow(dist, 2)) / 2);
    *result = exp(-1*power);
    return 0;
}

/**
 * returns the similiarity matrix A as defined in 1.1
 * @dataset_mat: the dataset X
 * @n: number of points in dataset
 * @dimension: dimension of points in the dataset
 * @sim_mat: out parameter. Will hold the similiarity matrix A.
 *  Return: 0 on success, 1 on failure
 */
int create_sim_mat(double **dataset_mat, int n, int dimension, double ***sim_mat){
    int i,j, failure;
    double aij_res;

    failure = allocate_double_matrix(n, n, sim_mat);
    CHECK_FAILURE(failure, error);

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j){ (*sim_mat)[i][j] = 0.0; }
            else{
                if (calc_aij_sim(dataset_mat[i], dataset_mat[j], dimension, &aij_res) != 0) { goto error; }
                (*sim_mat)[i][j] = aij_res;
            }
        }
    }
    return 0; /* OK! */

error:
    return 1;
}

/**
 * arr_sum - Calculates the sum of the elements in a double array.
 * @arr: the array
 * @n: the size of the array (needs to be predetermined)
 * @sum: outparameter. Will hold the sum of the array
 *  Return: 0 on success, 1 on failure
 */
int arr_sum(double* arr, int n, double* sum){
    int i;
    *sum = 0;

    for (i = 0; i < n; i++){
        *sum += arr[i];
    }

    return 0;
}

/**
 * create_diag_mat: Calculates the diagnoal degeree matrix as defined in 1.2
 * @dataset_mat: the input dataset X
 * @n: the length of dataset X
 * @dimension: dimension of points in dataset X
 * @diag_mat: out parameter. Will hold the diagonal degree matrix.
 *  Return: 0 on success, 1 on failure
 */
int create_diag_mat(double **dataset_mat, int n, int dimension, double ***diag_mat){
    int i, failure;
    double res;
    double** sim_mat = NULL;
    *diag_mat = NULL;

    if (create_sim_mat(dataset_mat, n, dimension, &sim_mat) != 0) { goto error; };

    failure = allocate_double_matrix(n,n, diag_mat);
    CHECK_FAILURE(failure, error);

    for (i = 0; i < n; i++){
        if (arr_sum(sim_mat[i], n, &res) != 0) { goto error; }
        (*diag_mat)[i][i] = res;
    }

    free_matrix(&sim_mat, n);
    sim_mat = NULL;
    return 0;  /* OK! */

error:
    if (sim_mat) { free_matrix(&sim_mat, n); }
    return 1;
}

/**
 * one_div_sqrt_diag_mat: Given a diagonal n*n matrix, Calculates the 1/sqrt(mat)
 * @mat: the matrix to calculate based on
 * @n: the dimension of the matrix
 * @ret_mat: out parameter. Will hold the result matrix
 *  Return: 0 on success, 1 on failure
 */
int one_div_sqrt_diag_mat(double** mat, int n, double ***ret_mat){
    int i, failure;
    *ret_mat = NULL; /* safe initialization of out parameter */
    
    failure = allocate_double_matrix(n,n, ret_mat);
    CHECK_FAILURE(failure, error);

    for (i = 0; i < n; i++){
        (*ret_mat)[i][i] = 1.0 / sqrt(mat[i][i]);
    }
    return 0;

error:
    return 1;
}

/**
 * diag_mat_mult - Calculates D*A matrix (direction == 'left') or A*D (direction == 'right')
 * @pre: D diagnoal n*n matrix
 * @D: dianoal n*n matrix
 * @A: n*n matrix
 * @direction: 'left' if D*A, 'right' if A*D
 * @res: out parameter. Will hold the result matrix
 *  Return: 0 on success, 1 on failure
*/
int diag_mat_mult(double** D, double** A, int n, char* direction, double*** res){
    int i,j, failure;
    *res = NULL; /* safe initialization for error handling */

    if(strcmp(direction, "left") != 0 && strcmp(direction, "right")) {goto error;}
    
    failure = allocate_double_matrix(n, n, res);
    CHECK_FAILURE(failure, error);

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (strcmp(direction, "left") == 0)
                (*res)[i][j] = D[i][i] * A[i][j];
            if (strcmp(direction, "right") == 0)
                (*res)[i][j] = A[i][j] * D[j][j];
        }
    }
    return 0;  /* OK! */
error:
    return 1;
}

/**
 * create_normalized_sim_mat - Calcualtes the graph Laplacian matrix, W, as defined in 1.3
 * @points_lst: the input dataset X
 * @n: the length of dataset X
 * @w_mat: outparameter. Will hold the result W matrix
 *  Return: 0 on success, 1 on failure
 */
int create_normalized_sim_mat(double **dataset_mat, int n, int dimension, double ***w_mat){
    double **sim_mat = NULL, **diag_mat_minus_sqrt = NULL, **diag_mat = NULL, **DA = NULL;

    if (create_sim_mat(dataset_mat, n, dimension, &sim_mat) != 0) { goto error; };
    if (create_diag_mat(dataset_mat, n, dimension, &diag_mat) != 0) { goto error; }
    if (one_div_sqrt_diag_mat(diag_mat, n, &diag_mat_minus_sqrt) != 0) { goto error; }

    if(diag_mat_mult(diag_mat_minus_sqrt, sim_mat, n, "left", &DA) != 0) { goto error; }
    if(diag_mat_mult(diag_mat_minus_sqrt, DA, n, "right", w_mat) != 0) { goto error; }

    /* free memory allocations */
    free_matrix(&sim_mat, n);
    free_matrix(&diag_mat, n);
    free_matrix(&diag_mat_minus_sqrt, n);
    free_matrix(&DA, n);
    return 0;

error:
    free_matrix(&sim_mat, n);
    free_matrix(&diag_mat, n);
    free_matrix(&diag_mat_minus_sqrt, n);
    free_matrix(&DA, n);
    return 1;
}

/**
 * frobenius_norm - Calculates the Frobenius norm as described in the Wikipedia page https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm
 * @mat: matrix of size n*m to calculate its Frobenius norm
 * @n: column dimension of mat
 * @m: row dimension of mat
 * @norm_res: out parameter. Will hold the Frobenius norm result.
 *  Return: 0 on success, 1 on failure
 */ 
int frobenius_norm(double **mat, int m, int n, double *norm_res){
    int i, j;
    *norm_res = 0;
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            *norm_res += (mat[i][j] * mat[i][j]);
        }
    }
    *norm_res = sqrt((*norm_res));
    return 0;
}

/**
 * substract_matrix - Given matrix A and B of size m*n each calculates A-B
 * @A: matrix A
 * @B: matrix B 
 * @m: number of rows
 * @n: number of columns
 * @res_mat: out parameter. Will hold the result matrix after the substraction
 *  Return: 0 on success, 1 on failure
 */
int substract_matrix(double **A, double **B, int m, int n, double ***res_mat){
    int i, j, failure;
    *res_mat = NULL;

    failure = allocate_double_matrix(m, n, res_mat);
    CHECK_FAILURE(failure, error);

    /* Preform matrix substraction */
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            (*res_mat)[i][j] = A[i][j] - B[i][j];
        }
    }

    return 0;

error:
    return 1;
}

/**
 * calc_h_next_ij - Calculates the next "iteration" of H matrix for for cell [i][j] as defined in 1.4.2
 * @i: the row index
 * @j: the column index
 * @wh: W*H matrix
 * @hh_th: H * H^t * H matrix
 * @res_ij: out parameter. Will hold the result of the calculation
 *  Return: 0 on success, 1 on failure
 */
int calc_h_next_ij(int i, int j, double **w_h, double **h_htr_h, double **h, double *res_ij){
    double beta = 0.5;  /* default value */
    (*res_ij) = h[i][j] * (1- beta + beta * (w_h[i][j] / h_htr_h[i][j]));
    return 0;
}

/**
 * multi_mat - multiplies 2 matrices, A (a*b) and B (b*c)
 * @a: row dimension of A
 * @b: column dimension of A + row dimension of B
 * @c: column dimension of B
 * @res_mat: out parameter. Will hold the result matrix of dimension a*c
 *  Return: 0 on success, 1 on failure
 */
int multi_mat(int a, int b, int c, double **A, double **B, double ***res_mat){
    int i, j, k, failure;
    double sum;
    *res_mat = NULL;

    failure = allocate_double_matrix(a, c, res_mat);
    CHECK_FAILURE(failure, error);

    /* Perform matrix multiplication */
    for (i = 0; i < a; i++){
        for (j = 0; j < c; j++){
            sum = 0;
            for (k = 0; k < b; k++){
                sum += A[i][k] * B[k][j];
            }
            (*res_mat)[i][j] = sum;
        }
    }
    return 0;

error:
    return 1;
}

/**
 * transpose_matrix - Transposes a matrix of size m x n into a matrix of size n x m
 * @mat: input matrix to transpose
 * @m: number of rows in input matrix
 * @n: number of columns in input matrix 
 * @transposed: out parameter. Will hold the transposed matrix
 * Return: 0 on success, 1 on failure
 */
int transpose_matrix(double **mat, int m, int n, double ***transposed){
    int i, j, failure;
    *transposed = NULL;

    failure = allocate_double_matrix(n, m, transposed); /* rows and cols dimensions swap becuase of transpose */
    CHECK_FAILURE(failure, error);

    /* preform the transpose operation */
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            (*transposed)[j][i] = mat[i][j];
        }
    }
    return 0;

error:
    return 1;
}

/**
 * calc_h_next - Calculates the next "iteration" of H matrix as defined in 1.4.2
 * @wh: W*H matrix
 * @hh_th: H * H^t * H matrix
 * @n: number of rows in the output matrix
 * @k: number of columns in the output matrix
 * @next_h: out parameter. Will hold the result of the calculation
 */
int calc_h_next(double **h, double **wh, double **hh_th, int n, int k, double ***next_h){
    int i, j, failure;
    double res;
    *next_h = NULL;

    failure = allocate_double_matrix(n, k, next_h);
    CHECK_FAILURE(failure, error);

    /* Perform the calcuations */
    for (i = 0; i < n; i++){
        for (j = 0; j < k; j++){
            if (calc_h_next_ij(i, j, wh, hh_th, h, &res) != 0) { goto error; }
            (*next_h)[i][j] = res;
        }
    }
    return 0;

error:
    if (*next_h) { free_matrix(next_h, n); }
    return 1;
}

/**
 * run the optimization algorithm as defined in 1.4
 * gets an input all the required matrices for calculations
 * @optimized_h: out parameter. Will hold the optimized H matrix
 */
int optimize_h_mat(double **h_init, double** w, int n, int k, double ***optimized_h){
    int i;
    double for_norm, **w_h, **htr, **h, **h_next, **h_htr, **h_htr_h, **h_next_minus_h_curr;
    h = h_init;
    for (i = 0; i < MAX_ITER; i++){
        CHECK_FAILURE(multi_mat(n, n, k, w, h, &w_h), error); /* W*H */
        CHECK_FAILURE(transpose_matrix(h, n, k, &htr), error);  /* H^t */
        CHECK_FAILURE(multi_mat(n, k, n, h, htr, &h_htr), error); /* H*H^t */
        CHECK_FAILURE(multi_mat(n, n, k, h_htr, h, &h_htr_h), error);/* H*H^t*H */
        CHECK_FAILURE(calc_h_next(h, w_h, h_htr_h, n, k, &h_next), error); /* H_(t+1) */
        CHECK_FAILURE(substract_matrix(h_next, h, n, k, &h_next_minus_h_curr), error); /* H_(t+1) - H_t */
        CHECK_FAILURE(frobenius_norm(h_next_minus_h_curr, n, k, &for_norm), error);
        free_matrix(&w_h, n);
        free_matrix(&htr, k);
        free_matrix(&h_htr, n);
        free_matrix(&h_htr_h, n);
        free_matrix(&h_next_minus_h_curr, n);
        if ((for_norm * for_norm) < EPSILON){
            (*optimized_h) = h_next;
            if (h != h_init) { free_matrix(&h, n); }
            return 0;  
        }
        if (h != h_init) free_matrix(&h,n); /* set h_t values to be h_next */
        h = h_next;
    }
    (*optimized_h) = h;
    return 0;
error:
    if (w_h) free_matrix(&w_h, n);
    if (htr) free_matrix(&htr, k);
    if (h_htr) free_matrix(&h_htr, n);
    if (h_htr_h) free_matrix(&h_htr_h, n);
    if (h_next_minus_h_curr) free_matrix(&h_next_minus_h_curr, n);
    return 1;
}

/**
 *  Validates command-line arguments and sets output parameters.
 * @argc:  Number of command-line arguments.
 * @argv:  Array of command-line argument strings.
 * @goal:  Out parameter. Will point to argv[1] if it is a valid goal ("sym", "norm", "ddg").
 * @path:  Out parameter. Will point to argv[2] if it is a valid existing file path.
*/
int validate_and_set_input(int argc, char *argv[], char **goal, char **path){
    FILE *file;
    const int EXPECTED_NUMBER_OF_ARGS = 3;
    *goal = NULL;
    *path = NULL;
    
    if (argc != EXPECTED_NUMBER_OF_ARGS) {goto error;}

    *goal = argv[1];
    *path = argv[2];

    /* check goal validity */
    if ((strcmp(*goal, "sym") != 0) && (strcmp(*goal, "norm") != 0) && (strcmp(*goal, "ddg") != 0)) {goto error;}

    /* check file exists by opening it for reading*/
    file = fopen(*path, "r");
    if (file == NULL) {goto error;}
    fclose(file);
    return 0;

error:
    return 1;
}

int main(int argc, char *argv[]){
    int failure, n, cols;
    char *goal, *path, *mat_str = NULL;
    double **ret_mat = NULL;
    double **dataset_mat = NULL;

    failure = validate_and_set_input(argc, argv, &goal, &path);
    CHECK_FAILURE(failure, error);

    failure = input_txt_to_matrix(path, &n, &cols, &dataset_mat);
    CHECK_FAILURE(failure, error);

    if (strcmp(goal, "sym") == 0){
        failure = create_sim_mat(dataset_mat, n, cols, &ret_mat);
        CHECK_FAILURE(failure, error);
    }
    else if (strcmp(goal, "ddg") == 0){
        failure = create_diag_mat(dataset_mat, n, cols, &ret_mat);
        CHECK_FAILURE(failure, error);
    }
    else if (strcmp(goal, "norm") == 0)
    {
        failure = create_normalized_sim_mat(dataset_mat, n, cols, &ret_mat);
        CHECK_FAILURE(failure, error);
    }

    failure = mat_to_str(ret_mat, n, n, &mat_str);
    CHECK_FAILURE(failure, error);

    printf("%s\n", mat_str);
    free(mat_str); /* free the matrix str*/
    free_matrix(&ret_mat, n);  /* free the matrix */
    free_matrix(&dataset_mat, n); /* free dataset matrix */
    return 0;
error:
    if (mat_str) free(mat_str);
    if (ret_mat) free_matrix(&ret_mat, n);
    if (dataset_mat) free_matrix(&dataset_mat, n);
    printf("%s", ERR_MSG_GENERAL);
    exit(1);
}