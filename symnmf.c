#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct point {
    struct point *next;
    struct cord *cords;
    int dim
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

    for (i = 0; i < dim1; i++){
        *dist += pow((curr1->value - curr2->value), 2);
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
    *result = exp(power);
    return 0;
}

/**
 * returns the similiarity matrix A as defined in 1.1
 * @points_lst: the dataset X
 * @sim_mat: out parameter. Will hold the similiarity matrix A.
 */
int create_sim_mat(point* points_lst, double*** sim_mat){
    int n = points_lst->dim;
    int i,j;
    double aij_res;

    point* head = points_lst;
    point* xi = points_lst;
    point* xj = points_lst;

    double** sim_mat = malloc(n * sizeof(double*));
    if (sim_mat == NULL) {return 1;}

    for (i = 0; i < n; i++){
        xj = head; 
        sim_mat[i] = malloc(n * sizeof(double));  /* allocate each row */
        if (sim_mat[i] == NULL) {return 1;}
        
        for (j = 0; j < n; j++){
            if (i == j){
                *sim_mat[i][j] = 0.0;
            }
            else{
                if (calc_aij_sim(xi, xj, &aij_res) != 0) {return 1;}
                *sim_mat[i][j] = aij_res;
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
 * @diag_mat: out parameter. Will hold the diagonal degree matrix.
 */
int create_diag_mat(point* points_lst, double ***diag_mat){
    int n = points_lst->dim;
    int i,j;
    double** sim_mat;
    double res;
    if (create_sim_mat(points_lst, &sim_mat) != 0) {return 1;};

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
    int tmp_val;
    *ret_mat = malloc(n * sizeof(double*));
    if (*ret_mat == NULL) {return 1;}
    for (i = 0; i < n; i++){
        (*ret_mat)[i] == calloc(n, sizeof(double));
        if ((*ret_mat[i]) == NULL) {return 1;}
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
    if (direction != "left" || direction != "right") {return 1;}

    int i,j;
    *res = malloc(n * sizeof(double*));
    if (*res == NULL) {return 1;}
    
    for (i = 0; i < n; i++){
        (*res)[i] = malloc(n * sizeof(n));
        if ((*res)[i] == NULL) {return 1;}

        for (j = 0; j < n; j++){
            if (direction == "left")
                (*res)[i][j] = D[i][i] * A[i][j];
            if (direction == "right")
                (*res)[i][j] = A[i][j] * D[j][j];
        }
    }
    return 0;  // OK!
}

/**
 * create_normalized_sim_mat - Calcualtes the graph Laplacian matrix, W, as defined in 1.3
 * @points_lst: the input dataset X
 * @w_mat: outparameter. Will hold the result W matrix
 */
int create_normalizec_sim_mat(point* points_lst, double*** w_mat){
    int n = points_lst->dim;
    double **sim_mat, **diag_mat_minus_sqrt, **diag_mat, **DA;

    if (create_sim_mat(points_lst, &sim_mat) != 0) {return 1;};
    if (create_diag_mat(points_lst, &diag_mat) != 0) {return 1;}
    if (one_div_sqrt_diag_mat(diag_mat, n, &diag_mat_minus_sqrt) != 0) {return 1;}

    if(diag_mat_mult(diag_mat_minus_sqrt, sim_mat, n, "left", &DA) != 0) {return 1;};
    if(diag_mat_mult(diag_mat_minus_sqrt, DA, n, "right", w_mat) != 0) {return 1;};

    return 0;
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