#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int calc_file_dimension(char *path, int *rows_out, int *cols_out){
    FILE *file;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int i, lines_count = 0, expected_cols = -1, curr_cols_count;

    file = fopen(path, "r");
    if (file == NULL) goto error;

    while ((read = getline(&line, &len, file)) != -1){
        if (read == 1 && line[0] == '\n') continue; /* last line is expected to be blank */

        lines_count++; /* increment number of lines in file */

        /* count columns in line */
        curr_cols_count = 1;
        for (i = 0; i < read; i++){
            if (line[i] == ',') curr_cols_count++;
        }

        /* if first line set the expected cols */
        if (expected_cols == -1) {
            expected_cols = curr_cols_count;
        } else if (curr_cols_count != expected_cols) goto error; /* inconsistent number of columns in file - error */
    }

    /* set out parameters */
    *rows_out = lines_count;
    *cols_out = expected_cols;
    
    fclose(file);
    return 0;

error:
    if (file) fclose(file);
    return 1;
}

int main(){
    int rows, cols;
    int failure;
    char *path = "./tests/input_1.txt";

    failure = calc_file_dimension(path, &rows, &cols);

    printf("failure = %d, rows = %d, cols = %d",failure, rows, cols);
    return 0;
}