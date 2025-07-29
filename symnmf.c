#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int failure;

typedef struct point {
    struct point *next;
    struct cord *cords;
    int dim;
} point;

typedef struct cord {
    double value;
    struct cord *next;
} cord;

//TODO
void general_error(){
    printf("An Error Has Occured\n");
    //Add freeing the memory
    exit(1);
}

/**
 * sqr_mat_to_str - returns a string describing the matrix
 * @mat: n*n matrix to be descibed
 * @n: the dimension of mat
 * @ret_str: outparameter. Will hold the describing string
 */
int sqr_mat_to_str(double** mat, int n, char **ret_str){
    int i, j, len = 0;
    int bufsize = n * n * 32 + n; // Estimate: 32 chars per number + newlines
    *ret_str = malloc(bufsize);
    if (*ret_str == NULL) {return 1;}

    (*ret_str)[0] = '\0'; // Start with empty string

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            // Append number, add comma if not last in row
            len += snprintf(*ret_str + len, bufsize - len, "%.4f", mat[i][j]);
            if (j < n - 1)
                len += snprintf((*ret_str) + len, bufsize - len, ",");
        }
        if (i < n - 1)
            len += snprintf((*ret_str) + len, bufsize - len, "\n");
    }
    return 0;
}

int input_txt_to_points_lst(char* path, point **head_point, int *points_count){
    double n;
    char c;
    int dim = 0;
    *points_count = 0;

    FILE *file = fopen(path, "r");  // Open for reading
    if (file == NULL) {return 1;}

    point *curr_point, *prev_point;
    cord *head_cord, *curr_cord;

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
        if (c == ','){ // still on same point
            dim++;
            curr_cord->value = n;
            curr_cord->next = (cord*) calloc(1, sizeof(cord));
            if (curr_cord->next == NULL) {return 1;}
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }

        if (c == '\n'){  // end of this point
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
        prev_point->next = NULL;
    }

    return 0;
}

/**
 * euc_distance - Calculates the eucledian distance between p1 and p2
 * @p1: First point
 * @p2: Second point
 * @dist: output parameter. Will hold the eucledian distance of the 2 points.
 */
int euc_distance(point *p1, point *p2, double *dist){
    int i;
    *dist = 0;
    if (p1->dim != p2->dim) {return 1;} // ERROR! two points not of same dimension

    int dim1 = p1->dim;
    cord *curr1, *curr2;

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
 */
int calc_aij_sim(point *xi, point *xj, double *result){
    *result = 0;
    double dist;
    if(euc_distance(xi, xj, &dist) != 0) {general_error();}
    double power = ((pow(dist, 2)) / 2);
    *result = exp(-1*power);
    return 0;
}

/**
 * returns the similiarity matrix A as defined in 1.1
 * @points_lst: the dataset X
 * @n: points_lst length
 * @sim_mat: out parameter. Will hold the similiarity matrix A.
 */
int create_sim_mat(point* points_lst, int n, double ***sim_mat){
    int i,j;
    double aij_res;

    point* head = points_lst;
    point* xi = points_lst;
    point* xj = points_lst;

    *sim_mat = malloc(n * sizeof(double*));
    if (*sim_mat == NULL) {return 1;}

    for (i = 0; i < n; i++){
        xj = head; 
        (*sim_mat)[i] = malloc(n * sizeof(double));  /* allocate each row */
        if ((*sim_mat)[i] == NULL) {return 1;}
        
        for (j = 0; j < n; j++){
            if (i == j){
                (*sim_mat)[i][j] = 0.0;
            }
            else{
                if (calc_aij_sim(xi, xj, &aij_res) != 0) {return 1;}
                (*sim_mat)[i][j] = aij_res;
            }
            xj = xj->next;
        }
        xi = xi->next;
    }

    return 0; // OK!
}

/**
 * arr_sum - Calculates the sum of the elements in a double array.
 * @arr: the array
 * @n: the size of the array (needs to be predetermined)
 * @sum: outparameter. Will hold the sum of the array
 */
int arr_sum(double* arr, int n, double* sum){
    *sum = 0;
    int i;

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
 */
int create_diag_mat(point* points_lst, int n, double ***diag_mat){
    int i,j;
    double** sim_mat;
    double res;
    if (create_sim_mat(points_lst, n, &sim_mat) != 0) {return 1;};

    *diag_mat = malloc(n * sizeof(double*));
    if (*diag_mat == NULL) {return 1;}

    for (i = 0; i < n; i++){
        (*diag_mat)[i] = calloc(n, sizeof(double));
        if ((*diag_mat)[i] == NULL) {return 1;}

        if (arr_sum(sim_mat[i], n, &res) != 0) {return 1;}
        (*diag_mat)[i][i] = res;
    }
    return 0;  // OK!
}

/**
 * one_div_sqrt_diag_mat: Given a diagonal n*n matrix, Calculates the 1/sqrt(mat)
 * @mat: the matrix to calculate based on
 * @n: the dimension of the matrix
 * @ret_mat: out parameter. Will hold the result matrix
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
*/
int diag_mat_mult(double** D, double** A, int n, char* direction, double*** res){
    if(strcmp(direction, "left") != 0 && strcmp(direction, "right")) {return 1;}

    int i,j;
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
    return 0;  // OK!
}

/**
 * create_normalized_sim_mat - Calcualtes the graph Laplacian matrix, W, as defined in 1.3
 * @points_lst: the input dataset X
 * @n: the length of dataset X
 * @w_mat: outparameter. Will hold the result W matrix
 */
int create_normalizec_sim_mat(point* points_lst, int n, double ***w_mat){
    double **sim_mat, **diag_mat_minus_sqrt, **diag_mat, **DA;

    if (create_sim_mat(points_lst, n, &sim_mat) != 0) {return 1;};
    if (create_diag_mat(points_lst, n, &diag_mat) != 0) {return 1;}
    if (one_div_sqrt_diag_mat(diag_mat, n, &diag_mat_minus_sqrt) != 0) {return 1;}

    if(diag_mat_mult(diag_mat_minus_sqrt, sim_mat, n, "left", &DA) != 0) {return 1;};
    if(diag_mat_mult(diag_mat_minus_sqrt, DA, n, "right", w_mat) != 0) {return 1;};

    return 0;
}

/**
 * frobenius_norm - Calculates the Frobenius norm as described in the Wikipedia page https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm
 * @mat: matrix of size n*m to calculate its Frobenius norm
 * @n: column dimension of mat
 * @m: row dimension of mat
 * @norm_res: out parameter. Will hold the Frobenius norm result.
 */ 
int frobenius_norm(double **mat, int n, int m, double *norm_res){
    *norm_res = 0;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            *norm_res += pow(mat[i][j], 2);
        }
    }
    *norm_res = sqrt((*norm_res));
    return 0;
}

/**
 * avg_mat_val - Calculates a matrice's avg value of all cells
 * @mat: n*m matrix (cells are of type double)
 * @n: columns dimension
 * @m: rows dimension
 * @avg_val: out parameter. Will hold the average value of mat's cells.
 */
int avg_mat_val(double **mat, int n, int m, double *avg_val){
    *avg_val = 0;

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            (*avg_val) += mat[i][j];
        }
    }

    return 0;  // OK!
}

/**
 * calc_h_next_ij - Calculates the next "iteration" of H matrix for for cell [i][j] as defined in 1.4.2
 * @i: the row index
 * @j: the column index
 * @wh: W*H matrix
 * @hh_th: H * H^t * H matrix
 * @res_ij: out parameter. Will hold the result of the calculation
 */
int calc_h_next_ij(int i, int j, double **wh, double **hh_th, double *res_ij){
    double beta = 0.5;  // Default value
    double mat_calc = (wh[i][j]) / (hh_th[i][j]);
    (*res_ij) = 1 - beta + (beta * mat_calc);

    return 0;
}

/**
 * multi_mat - multiplies 2 matrices, A (a*b) and B (b*c)
 * @a: row dimension of A
 * @b: column dimension of A + row dimension of B
 * @c: column dimension of B
 * @res_mat: out parameter. Will hold the result matrix of dimension a*c
 */
int multi_mat(int a, int b, int c, double **A, double **B, double ***res_mat){
    int i,j,k;
    double sum;

    // Allocate memory for result matrix
    *res_mat = malloc(a * sizeof(double));
    if (*res_mat == NULL) {return 1;}

    for (i = 0; i < a; i++){
        (*res_mat)[i] = calloc(c, sizeof(double));
        if ((*res_mat)[i] == NULL) {return 1;}
    }

    // Perform matrix multiplication
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
 * calc_h_next - Calculates the next "iteration" of H matrix as defined in 1.4.2
 * @wh: W*H matrix
 * @hh_th: H * H^t * H matrix
 * @n: number of rows in the output matrix
 * @k: number of columns in the output matrix
 * @next_h: out parameter. Will hold the result of the calculation
 */
int calc_h_next(double **wh, double **hh_th, int n, int k, double ***next_h){
    double res;
    (*next_h) = malloc(n * sizeof(double*));
    if (*next_h == NULL) {return 1;}

    // Allocate memory for result matrix
    for (int i = 0; i < n; i++){
        (*next_h)[i] = calloc(k, sizeof(double));
        if ((*next_h)[i] == NULL) {return 1;}
    }

    // Perform the calcuations
    for (int i = 0; i < n; i++){
        for (int j = 0; j < k; j++){
            if (calc_h_next_ij(i, j, wh, hh_th, &res) != 0) {return 1;}
            (*next_h)[i][j] = res;
        }
    }

    return 0;
}
