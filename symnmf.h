#ifndef SYMNMF_H
#define SYMNMF_H

# include <string.h>

typedef struct point {
    struct point *next;
    struct cord *cords;
    int dim;
} point;

typedef struct cord {
    double value;
    struct cord *next;
} cord;

int free_cords(cord *crd);

int free_pnt_lst(point *pnt);

int free_matrix(double **mat, int m);

int create_sim_mat(point* points_lst, int n, double ***sim_mat);

int mat_to_str(double** mat, int m, int n, char **ret_str);

int create_diag_mat(point* points_lst, int n, double ***diag_mat);

int create_normalized_sim_mat(point* points_lst, int n, double ***w_mat);

int optimize_h_mat(double **init_h_mat, double** w, int n, int k, double ***optimized_h);

#endif