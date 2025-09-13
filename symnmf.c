#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CHECK_FAILURE(failure, label) \
    do { if ((failure) == 1) goto label; } while (0)

char *ERR_MSG_GENERAL = "An Error Has Occurred\n";
int failure;
double EPSILON = 0.0001;
int MAX_ITER = 300;

typedef struct point {
    struct point *next;
    struct cord *cords;
    int dim;
} point;

typedef struct cord {
    double value;
    struct cord *next;
} cord;


/**
 * free_cords - free cordinates memory allocation for the linked list starting at crd
 * @crd: the starting node to free from
 *  Return: 0 on success, 1 on failure
 */
int free_cords(cord *crd){
    cord *nxt;

    while (crd){
        nxt = crd->next;
        free(crd);
        crd = nxt;
    }

    return 0;
}

/**
 * free_pnt_lst - free memory allocation for all points in the linked list (starting from pnt)
 * @pnt: the head node to start from
 *  Return: 0 on success, 1 on failure
*/ 
int free_pnt_lst(point *pnt){
    point *nxt;

    while (pnt){
        if (free_cords(pnt->cords) != 0) {return 1;}
        nxt = pnt->next;
        free(pnt);
        pnt = nxt;
    }

    return 0;
}

/**
 * free_matrix - Frees the memory allocated for a 2d matrix mat with m rows
 * @mat: the matrix (a 2d array)
 * @m: number of rows
 *  Return: 0 on success, 1 on failure
 */
int free_matrix(double **mat, int m){
    int i;

    if (mat == NULL) return 0;  /* nothing to free */

    /* Free each row */
    for (i = 0; i < m; i++){
        if (mat[i] != NULL)
            free(mat[i]);
    }

    free(mat);
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
 * input_txt_to_points_lst - Creates a linked list of points (and cords) for the .txt dataset given as an input
 * @path: the relative path to the input file
 * @head_point: out parameter. Will hold the first node (of type point) in the linked list.
 * @points_count: the number of points in the dataset = the dataset length
 *  Return: 0 on success, 1 on failure
 */
int input_txt_to_points_lst(char* path, point **head_point, int *points_count){
    double n;
    char c;
    FILE *file;
    point *curr_point, *prev_point = NULL;
    cord *head_cord, *curr_cord;
    int dim = 0;
    *points_count = 0;

    file = fopen(path, "r");  /* Open for reading */
    if (file == NULL) {return 1;}

    /* init head of cords*/
    head_cord = (cord*) malloc(sizeof(struct cord));
    if (head_cord == NULL) {return 1;}

    curr_cord = head_cord;
    curr_cord->next = NULL;

    /* init head of points linked_list */
    (*head_point) = (point*) malloc(sizeof(point));
    if (*head_point == NULL) {return 1;}

    curr_point = (*head_point);
    curr_point->next = NULL;


    /* scan input */
    while (fscanf(file, "%lf%c", &n, &c) == 2){
        if (c == ','){ /* still on same point */
            dim++;
            curr_cord->value = n;
            curr_cord->next = (cord*) calloc(1, sizeof(cord));
            if (curr_cord->next == NULL) {return 1;}
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }

        if (c == '\n'){  /* end of this point */
            curr_cord->value = n;
            curr_point->cords = head_cord;
            curr_point->dim = dim + 1;
            prev_point = curr_point;
            curr_point->next = (point*) calloc(1, sizeof(point));
            if (curr_point->next == NULL) {return 1;}

            curr_point = curr_point->next;
            curr_point->next = NULL;
            curr_point->cords = NULL;

            head_cord = (cord*) calloc(1, sizeof(cord));
            if (head_cord == NULL) {return 1;}
            curr_cord = head_cord;
            curr_cord->next = NULL;
            dim = 0;

            (*points_count)++;
        }
    }

    /* freeing memory for last point and its cord */
    if (curr_point != NULL){
        free(curr_point);
        free(head_cord);
        if (prev_point != NULL){
            prev_point->next = NULL;
        }
    }

    fclose(file);

    return 0;
}

/**
 * euc_distance - Calculates the eucledian distance between p1 and p2
 * @p1: First point
 * @p2: Second point
 * @dist: output parameter. Will hold the eucledian distance of the 2 points.
 *  Return: 0 on success, 1 on failure
 */
int euc_distance(point *p1, point *p2, double *dist){
    int i, dim1;
    cord *curr1, *curr2;
    *dist = 0;
    if (p1->dim != p2->dim) {return 1;} /* ERROR! two points not of same dimension */

    dim1 = p1->dim;
    curr1 = p1->cords;
    curr2 = p2->cords;

    for (i = 0; i < dim1; i++){
        (*dist) += pow((curr1->value - curr2->value), 2);
        curr1 = curr1->next;
        curr2 = curr2->next;
    }

    *dist = sqrt(*dist);
    return 0;
}

/**
 * calc_aij_sim - Calculates the value of A(ij) in the similiarity matrix as defined in 1.1
 * @pre: xi.dim == x2.dim
 * @xi: the first cordinate in the dataset X
 * @xj: the second cordiante in the dataset X 
 * @result: output parameter. Will hold the result of A(ij)
 *  Return: 0 on success, 1 on failure
 */
int calc_aij_sim(point *xi, point *xj, double *result){
    double dist, power;
    *result = 0;
    if(euc_distance(xi, xj, &dist) != 0) { return 1;}
    power = ((pow(dist, 2)) / 2);
    *result = exp(-1*power);
    return 0;
}

/**
 * returns the similiarity matrix A as defined in 1.1
 * @points_lst: the dataset X
 * @n: points_lst length
 * @sim_mat: out parameter. Will hold the similiarity matrix A.
 *  Return: 0 on success, 1 on failure
 */
int create_sim_mat(point* points_lst, int n, double ***sim_mat){
    int i,j,k;
    double aij_res;

    point* head = points_lst;
    point* xi = points_lst;
    point* xj = points_lst;

    *sim_mat = malloc(n * sizeof(double*));
    if (*sim_mat == NULL) {return 1;}

    /* allocate rows */ 
    for (i = 0; i < n; i++){
        xj = head; 
        (*sim_mat)[i] = malloc(n * sizeof(double));  /* allocate each row */
        if ((*sim_mat)[i] == NULL) {
            /* free previousley allocated rows */
            for (j = 0; j < i; j++){
                free((*sim_mat)[j]);
            }
            free(*sim_mat);
            return 1;
        }
        
        /* populate the matrix */
        for (j = 0; j < n; j++){
            if (i == j){
                (*sim_mat)[i][j] = 0.0;
            }
            else{
                if (calc_aij_sim(xi, xj, &aij_res) != 0) { /* Error */
                    /* free all allocated rows*/
                    for (k = 0; k <= i; k++){ free((*sim_mat)[k]); }
                    free(*sim_mat);
                    return 1;
                }
                (*sim_mat)[i][j] = aij_res;
            }
            xj = xj->next;
        }
        xi = xi->next;
    }

    return 0; /* OK! */

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
 * @points_lst: the input dataset X
 * @n: the length of dataset X
 * @diag_mat: out parameter. Will hold the diagonal degree matrix.
 *  Return: 0 on success, 1 on failure
 */
int create_diag_mat(point* points_lst, int n, double ***diag_mat){
    int i,j;
    double** sim_mat;
    double res;
    if (create_sim_mat(points_lst, n, &sim_mat) != 0) {return 1;};

    *diag_mat = malloc(n * sizeof(double*));
    if (*diag_mat == NULL) {
        /* free allocated memory */
        free_matrix(sim_mat, n);
        return 1;
    }

    for (i = 0; i < n; i++){
        (*diag_mat)[i] = calloc(n, sizeof(double));
        if ((*diag_mat)[i] == NULL) {
            /* free previousley allocated memory */
            for (j = 0; j < i; j++){ free((*diag_mat)[j]); }
            free(*diag_mat);
            free_matrix(sim_mat, n);
            return 1;
        }

        if (arr_sum(sim_mat[i], n, &res) != 0) {
            /* free previosley allocated rows */
            for (j = 0; j < i; j++) { free((*diag_mat)[j]); }
            free(*diag_mat);
            free_matrix(sim_mat, n);
            return 1;
        }
        
        (*diag_mat)[i][i] = res;
    }

    free_matrix(sim_mat, n);
    return 0;  /* OK! */
}

/**
 * one_div_sqrt_diag_mat: Given a diagonal n*n matrix, Calculates the 1/sqrt(mat)
 * @mat: the matrix to calculate based on
 * @n: the dimension of the matrix
 * @ret_mat: out parameter. Will hold the result matrix
 *  Return: 0 on success, 1 on failure
 */
int one_div_sqrt_diag_mat(double** mat, int n, double ***ret_mat){
    int i;
    (*ret_mat) = malloc(n * sizeof(double*));
    if ((*ret_mat) == NULL) {return 1;}
    for (i = 0; i < n; i++){
        (*ret_mat)[i] = calloc(n, sizeof(double));
        if ((*ret_mat)[i] == NULL) {return 1;}
        (*ret_mat)[i][i] = 1.0 / sqrt(mat[i][i]);
    }
    return 0;
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
    int i,j;
    if(strcmp(direction, "left") != 0 && strcmp(direction, "right")) {return 1;}

    (*res) = malloc(n * sizeof(double*));
    if (*res == NULL) {return 1;}
    
    for (i = 0; i < n; i++){
        (*res)[i] = malloc(n * sizeof(double));
        if ((*res)[i] == NULL) {return 1;}

        for (j = 0; j < n; j++){
            if (strcmp(direction, "left") == 0)
                (*res)[i][j] = D[i][i] * A[i][j];
            if (strcmp(direction, "right") == 0)
                (*res)[i][j] = A[i][j] * D[j][j];
        }
    }
    return 0;  /* OK! */
}

/**
 * create_normalized_sim_mat - Calcualtes the graph Laplacian matrix, W, as defined in 1.3
 * @points_lst: the input dataset X
 * @n: the length of dataset X
 * @w_mat: outparameter. Will hold the result W matrix
 *  Return: 0 on success, 1 on failure
 */
int create_normalized_sim_mat(point* points_lst, int n, double ***w_mat){
    double **sim_mat, **diag_mat_minus_sqrt, **diag_mat, **DA;

    if (create_sim_mat(points_lst, n, &sim_mat) != 0) {return 1;};
    if (create_diag_mat(points_lst, n, &diag_mat) != 0) {return 1;}
    if (one_div_sqrt_diag_mat(diag_mat, n, &diag_mat_minus_sqrt) != 0) {return 1;}

    if(diag_mat_mult(diag_mat_minus_sqrt, sim_mat, n, "left", &DA) != 0) {return 1;};
    if(diag_mat_mult(diag_mat_minus_sqrt, DA, n, "right", w_mat) != 0) {return 1;};

    /* free memory allocations */
    free_matrix(sim_mat, n);
    free_matrix(diag_mat, n);
    free_matrix(diag_mat_minus_sqrt, n);
    free_matrix(DA, n);
    return 0;
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
    int i, j;

    /* Allocate memory for result matrix */
    *res_mat = malloc(m * sizeof(double*));
    if (*res_mat == NULL) return 1;

    /* Allocate memory for each row */
    for (i = 0; i < m; i++){
        (*res_mat)[i] = malloc(n * sizeof(double));
        if ((*res_mat)[i] == NULL){
            /* clean up already allocated memory */
            for (j = 0; j < i; j++){
                free((*res_mat)[j]);
            }
            free(*res_mat);
            return 1;
        }
    }

    /* Preform matrix substraction */
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            (*res_mat)[i][j] = A[i][j] - B[i][j];
        }
    }

    return 0;
}

/**
 * avg_mat_val - Calculates a matrice's avg value of all cells
 * @mat: n*m matrix (cells are of type double)
 * @n: columns dimension
 * @m: rows dimension
 * @avg_val: out parameter. Will hold the average value of mat's cells.
 *  Return: 0 on success, 1 on failure
 */
int avg_mat_val(double **mat, int n, int m, double *avg_val){
    int i, j;
    *avg_val = 0;

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            (*avg_val) += mat[i][j];
        }
    }

    return 0;  /* OK! */
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
    int i,j,k;
    double sum;

    /* Allocate memory for result matrix */
    *res_mat = malloc(a * sizeof(double));
    if (*res_mat == NULL) {return 1;}

    for (i = 0; i < a; i++){
        (*res_mat)[i] = calloc(c, sizeof(double));
        if ((*res_mat)[i] == NULL) {return 1;}
    }

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
    int i, j;

    /* Allocate memory for transposed matrix (n rows since dimensions are swapped) */
    *transposed = (double**) malloc(n * sizeof(double*));
    if (*transposed == NULL) return 1;

    /* Allocate memory for each row (m columns since dimesnions are swapped) */
    for (i = 0; i < n; i++){
        (*transposed)[i] = malloc(m * sizeof(double));
        if ((*transposed)[i] == NULL){
            /* clean up already allocated memory */
            for (j = 0; j < i; j++){
                free((*transposed)[j]);
            }
            free(*transposed);
            return 1;
        }
    }

    /* preform the transpose operation */
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            (*transposed)[j][i] = mat[i][j];
        }
    }

    return 0;
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
    int i, j;
    double res;
    (*next_h) = malloc(n * sizeof(double*));
    if (*next_h == NULL) {return 1;}

    /* Allocate memory for result matrix */
    for (i = 0; i < n; i++){
        (*next_h)[i] = calloc(k, sizeof(double));
        if ((*next_h)[i] == NULL) {return 1;}
    }

    /* Perform the calcuations */
    for (i = 0; i < n; i++){
        for (j = 0; j < k; j++){
            if (calc_h_next_ij(i, j, wh, hh_th, h, &res) != 0) {return 1;}
            (*next_h)[i][j] = res;
        }
    }

    return 0;
}

/**
 * run the optimization algorithm as defined in 1.4
 * @h_init: the initial H matrix
 * @w: W matrix
 * @n: size of dataset (= number of rows in H matrix)
 * @k: number of clusters (= number of columns in H matrix)
 * @optimized_h: out parameter. Will hold the optimized H matrix
 *  Return: 0 on success, 1 on failure
 */
int optimize_h_mat(double **h_init, double** w, int n, int k, double ***optimized_h){
    int i, failure;
    double for_norm;
    double **w_h, **htr, **h, **h_next, **h_htr, **h_htr_h, **h_next_minus_h_curr;
    
    h = h_init;
    for (i = 0; i < MAX_ITER; i++){
        /* calculte the needed matrices */
        failure = multi_mat(n, n, k, w, h, &w_h);  /* W*H */
        CHECK_FAILURE(failure, error);
        failure = transpose_matrix(h, n, k, &htr);  /* H^t */
        CHECK_FAILURE(failure, error);
        failure = multi_mat(n, k, n, h, htr, &h_htr); /* H*H^t */
        CHECK_FAILURE(failure, error);
        failure = multi_mat(n, n, k, h_htr, h, &h_htr_h);  /* H*H^t*H */
        CHECK_FAILURE(failure, error);
        failure = calc_h_next(h, w_h, h_htr_h, n, k, &h_next);  /* H_(t+1) */
        CHECK_FAILURE(failure, error);
        failure = substract_matrix(h_next, h, n, k, &h_next_minus_h_curr); /* H_(t+1) - H_t */
        CHECK_FAILURE(failure, error);
        failure = frobenius_norm(h_next_minus_h_curr, n, k, &for_norm);
        CHECK_FAILURE(failure, error);

        free_matrix(w_h, n);
        free_matrix(htr, k);
        free_matrix(h_htr, n);
        free_matrix(h_htr_h, n);
        free_matrix(h_next_minus_h_curr, n);

        if ((for_norm * for_norm) < EPSILON){
            (*optimized_h) = h_next;
            if (h != h_init) { free_matrix(h, n); }
            return 0;  
        }

        /* set h_t values to be h_next */
        if (h != h_init) free_matrix(h,n);
        h = h_next;
        CHECK_FAILURE(failure, error);
    }

    (*optimized_h) = h;
    return 0;

error:
    if (w_h) free_matrix(w_h, n);
    if (htr) free_matrix(htr, k);
    if (h_htr) free_matrix(h_htr, n);
    if (h_htr_h) free_matrix(h_htr_h, n);
    if (h_next_minus_h_curr) free_matrix(h_next_minus_h_curr, n);
    return 1;
}

int main(int argc, char *argv[]){
    int failure, n;
    char *goal, *file_name, *mat_str;
    double **ret_mat;
    point *dataset;

    if (argc != 3){
        printf("%s", ERR_MSG_GENERAL);
        return 1;
    }

    goal = argv[1];
    file_name = argv[2];

    failure = input_txt_to_points_lst(file_name, &dataset, &n);
    CHECK_FAILURE(failure, error);

    if (strcmp(goal, "sym") == 0){
        failure = create_sim_mat(dataset, n, &ret_mat);
        CHECK_FAILURE(failure, error);
    }
    else if (strcmp(goal, "ddg") == 0){
        failure = create_diag_mat(dataset, n, &ret_mat);
        CHECK_FAILURE(failure, error);
    }
    else if (strcmp(goal, "norm") == 0)
    {
        failure = create_normalized_sim_mat(dataset, n, &ret_mat);
        CHECK_FAILURE(failure, error);
    }
    else{
        printf("%s", ERR_MSG_GENERAL);
        goto error;
    }

    failure = mat_to_str(ret_mat, n, n, &mat_str);
    CHECK_FAILURE(failure, error);

    printf("%s\n", mat_str);

    /* free memory allocations */
    free(mat_str); /* free the matrix str*/
    free_matrix(ret_mat, n);  /* free the matrix */
    free_pnt_lst(dataset);  /* free the data set points */
    return 0;

error:
    if (mat_str) free(mat_str);
    if (ret_mat) free_matrix(ret_mat, n);
    if (dataset) free_pnt_lst(dataset);
    printf("%s", ERR_MSG_GENERAL);
    return 1;
}