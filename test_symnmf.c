#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

/**
 * void test_sym_goal(){
    point *dataset;
    int n, failure;

    failure = input_txt_to_points_lst(path_short, &dataset, &n);
    if (failure != 0) {printf("result is %d\n", failure);}

    double **sim_mat;
    failure = create_sim_mat(dataset, n, &sim_mat);
    if (failure != 0) {printf("result is %d\n", failure);}

    char* res;
    failure = mat_to_str(sim_mat, n, n, &res);
    if (failure != 0) {printf("result is %d\n", failure);}

    printf("%s\n", res);
}

void test_ddg_goal(){
    point *dataset;
    int n, failure;

    failure = input_txt_to_points_lst(path_short, &dataset, &n);
    if (failure != 0) {printf("result is %d\n", failure);}

    double **ddg_mat;
    failure = create_diag_mat(dataset, n, &ddg_mat);
    
    char* res;
    failure = mat_to_str(ddg_mat, n, n, &res);
    if (failure != 0) {printf("result is %d\n", failure);}

    printf("%s\n", res);
}

void test_normalized_similiarity_matrix(){
    point *dataset;
    int n, failure;

    failure = input_txt_to_points_lst(path_short, &dataset, &n);
    if (failure != 0) {printf("result is %d\n", failure);}

    double **norm_mat;
    failure = create_normalizec_sim_mat(dataset, n, &norm_mat);
    
    char* res;
    failure = mat_to_str(norm_mat, n, &res);
    if (failure != 0) {printf("result is %d\n", failure);}

    printf("%s\n", res);
}
 */

void test_file_parsing(){
    double **ret_mat;
    char *path = "./tests/input_1.txt";
    int rows, cols, failure;
    char *mat_str;

    failure = input_txt_to_matrix(path, &rows, &cols, &ret_mat);
    if (failure != 0) {printf("result is %d\n", failure);}

    failure = mat_to_str(ret_mat, rows, cols, &mat_str);
    if (failure != 0) {printf("result is %d\n", failure);}

    printf("rows = %d, cols = %d, and the matrix is:\n%s\n",rows, cols, mat_str);
}

int main() {
    test_file_parsing();
}