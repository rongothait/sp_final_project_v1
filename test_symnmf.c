#include <stdio.h>
#include <stdlib.h>
#include "symnmf.c" // Or use a header if you have one

char *path = "./input_1.txt";
char *path_short = "./input_1_short.txt";

void test_sym_goal(){
    point *dataset;
    int n, failure;

    failure = input_txt_to_points_lst(path_short, &dataset, &n);
    if (failure != 0) {printf("result is %d\n", failure);}

    double **sim_mat;
    failure = create_sim_mat(dataset, n, &sim_mat);
    if (failure != 0) {printf("result is %d\n", failure);}

    char* res;
    failure = sqr_mat_to_str(sim_mat, n, &res);
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
    failure = sqr_mat_to_str(ddg_mat, n, &res);
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
    failure = sqr_mat_to_str(norm_mat, n, &res);
    if (failure != 0) {printf("result is %d\n", failure);}

    printf("%s\n", res);
}

int main() {
    test_ddg_goal();
    printf("%c",'\n');
    
    test_sym_goal();
    printf("%c",'\n');

    test_normalized_similiarity_matrix();
}