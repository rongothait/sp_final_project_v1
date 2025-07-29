#include <stdio.h>
#include <stdlib.h>
#include "symnmf.c" // Or use a header if you have one

int main() {
    int n = 3;
    int i, j;
    double **mat;
    char *mat_str;

    // Allocate and fill a 3x3 matrix
    mat = malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        mat[i] = malloc(n * sizeof(double));
        for (j = 0; j < n; j++) {
            mat[i][j] = i * n + j + 1; // Fill with 1..9
        }
    }

    // Test sqr_mat_to_str
    if (sqr_mat_to_str(mat, n, &mat_str) != 0) {
        printf("Error in sqr_mat_to_str\n");
        return 1;
    }
    printf("Matrix as string:\n%s\n", mat_str);

    // Free memory
    free(mat_str);
    for (i = 0; i < n; i++) free(mat[i]);
    free(mat);

    return 0;
}