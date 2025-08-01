#ifndef SYMNMF_H
#define SYMNMF_H

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

int create_sim_mat(point* points_lst, int n, double ***sim_mat);

int sqr_mat_to_str(double** mat, int n, char **ret_str);

int create_diag_mat(point* points_lst, int n, double ***diag_mat);

#endif