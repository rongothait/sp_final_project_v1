#ifndef SYMNMF_H
#define SYMNMF_H

# include <string.h>

int free_matrix(double ***pmat, int m);

int create_sim_mat(double **dataset_mat, int n, int dimension, double ***sim_mat);

int mat_to_str(double** mat, int m, int n, char **ret_str);

int create_diag_mat(double **dataset_mat, int n, int dimension, double ***diag_mat);

int create_normalized_sim_mat(double **dataset_mat, int n, int dimension, double ***w_mat);

int optimize_h_mat(double **init_h_mat, double** w, int n, int k, double ***optimized_h);

int allocate_double_matrix(int rows, int cols, double ***ret_mat);

int calc_file_dimension(char *path, int *rows_out, int *cols_out);

int input_txt_to_matrix(char *path, int *rows_out, int *cols_out, double ***dataset_mat);

#endif