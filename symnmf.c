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

/* Calculates the eucledian distance between p1 and p2
   @pre: p1.dim = p2.dim
*/
double euc_distance(point *p1, point *p2){
    int i;
    double dist = 0;
    int dim1 = p1->dim;
    cord *curr1, *curr2;

    for (i = 0; i < dim1; i++){
        dist += pow((curr1->value - curr2->value), 2);

        curr1 = curr1->next;
        curr2 = curr2->next;
    }

    return sqrt(dist);
}

/*
calculates the value of A(ij) in the similarity matrix as defined
@pre: xi.dim = xj.dim
*/
double calc_aij_sim(point *xi, point *xj){
    double norm_sq = pow(euc_distance(xi, xj), 2);
    double power = - (norm_sq / 2);
    return exp(power);
}


double** create_sim_mat(point* points_lst){
    int n = points_lst->dim;
    int i,j;

    point* head = points_lst;
    point* xi = points_lst;
    point* xj = points_lst;

    double** sim_mat = malloc(n * sizeof(double*));
    if (sim_mat == NULL)
        general_error();

    for (i = 0; i < n; i++){
        xj = head; 
        sim_mat[i] = malloc(n * sizeof(double));  /* allocate each row */
        if (sim_mat[i] == NULL)
            general_error();
        
        for (j = 0; j < n; j++){
            if (i == j){
                sim_mat[i][j] == 0.0;
            }
            else{
                sim_mat[i][j] = calc_aij_sim(xi, xj);
            }
            xj = xj->next;
        }
        xi = xi->next;
    }

    return sim_mat;
}

/*
returns the sum of the array 'arr' of size 'n'
*/
double arr_sum(double* arr, int n){
    double sum = 0;
    int i;

    for (i = 0; i < n; i++){
        sum += arr[i];
    }

    return sum;
}

/*
returns the diagonal degree matrix
*/
double ** create_diag_mat(point* points_lst){
    int n = points_lst->dim;
    int i,j;

    double** sim_mat = create_sim_mat(points_lst);

    double** diag_mat = malloc(n * sizeof(double*));
    if (diag_mat == NULL)
        general_error();

    for (i = 0; i < n; i++){
        diag_mat[i] = calloc(n, sizeof(double));
        if (diag_mat[i] == NULL)
            general_error();

        diag_mat[i][i] = arr_sum(sim_mat[i], n);
    }

    return diag_mat;
}

/*
given a diagnoal n*n matrix, return 1/sqrt(mat)
*/
double** one_div_sqrt_diag_mat(double** mat, int n){
    int i;
    int tmp_val;
    double** ret_mat = malloc(n * sizeof(double*));
    if (ret_mat == NULL)
        general_error();
    for (i = 0; i < n; i++){
        ret_mat[i] == calloc(n, sizeof(double));
        if (ret_mat[i] == NULL)
            general_error();
        ret_mat[i][i] = 1.0 / sqrt(mat[i][i]);
    }

    return ret_mat;
}

/*
pre: D diagonal n*n matrix
@pre: A n*n matrix
returns the result matrix D*A (direction = 'left') or D*A (direction = 'right')
*/
double** diag_mat_mult(double** D, double** A, int n, char* direction){
    if (direction != "left" || direction != "right")
        general_error();

    int i,j;
    double **res = malloc(n * sizeof(double*));
    if (res == NULL)
        general_error();
    
    for (i = 0; i < n; i++){
        res[i] = malloc(n * sizeof(n));
        if (res[i] == NULL)
            general_error();

        for (j = 0; j < n; j++){
            if (direction == "left")
                res[i][j] = D[i][i] * A[i][j];
            if (direction == "right")
                res[i][j] = A[i][j] * D[j][j];
        }
    }

    return res;
}


double** create_nor_sim_mat(point* points_lst){
    int n = points_lst->dim;
    double** sim_mat = create_sim_mat(points_lst);
    double** diag_mat_minus_sqrt = one_div_sqrt_diag_mat(create_diag_mat(points_lst), n);

    double** DA = diag_mat_mult(diag_mat_minus_sqrt, sim_mat, n, "left");
    double** DAD = diag_mat_mult(diag_mat_minus_sqrt, DA, n, "right");

    return DAD;
}


char* sqr_mat_to_str(double** mat, int n){
    int i, j, len = 0;
    int bufsize = n * n * 32 + n; // Estimate: 32 chars per number + newlines
    char* buf = malloc(bufsize);
    if (buf == NULL)
        general_error();

    buf[0] = '\0'; // Start with empty string

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            // Append number, add comma if not last in row
            len += snprintf(buf + len, bufsize - len, "%.4f", mat[i][j]);
            if (j < n - 1)
                len += snprintf(buf + len, bufsize - len, ",");
        }
        if (i < n - 1)
            len += snprintf(buf + len, bufsize - len, "\n");
    }
    return buf;
}