import numpy as np
import sys
import symnmf_mod as symnmf
import numpy as np

ERR_MSG = "An Error Has Occured"

# set the seed as requested
np.random.seed(1234)

def general_error(msg = ERR_MSG):
    print(msg)
    exit(1)

"""
is_float - Checks if a string x represents a floating number
@x: the string
"""
def is_float(x):
    try:
        float(x)
        return True
    except:
        return False

"""
get_average_val_of_matrix - Calculates the average value of all of the cells in matrix
@mat: the matrix
"""
def get_average_val_of_matrix(mat):
    sum = 0
    for row in mat:
        for val in row:
            sum += val
    cnt = len(mat) * len(mat[0])
    return sum / cnt

"""
set_data - checks cmd args valididty
@data: an array of cmd args
"""
def set_data(data):
    # 1. making sure the input is in correct length
    if len(data) != 4:
        general_error()

    # 2. read arguments from CMD
    k = int(sys.argv[1])  # int, number of clusters
    goal = sys.argv[2]  # 'symnf' / 'sym' / 'ddg' / 'norm'
    path = sys.argv[3]  # path to data set


    # 3. check validity
    if not is_float(k):
        general_error()
    
    k = float(k)
    if (k <= 1 or k % 1 != 0):
        general_error()

    return int(k), goal, path

"""
txt_input_to_list - converts the path given to a list of lists
@path: the path of the file
"""
def txt_input_to_list(path):
    try:
        with open(path, 'r') as f:
            lines = [line.strip() for line in f]
            lines = [[float(x) for x in line.split(',')] for line in lines]
    except:
        general_error()
    
    return lines

"""
matrix_to_str - converts a list of lists (float) to formatted str ready for printing
@mat: the matrix to convert
"""
def matrix_to_str(mat):
    rows = []
    for row in mat:
        row_str = ','.join([f'{x:.4f}' for x in row])
        rows.append(row_str)

    return '\n'.join(rows)

def create_random_matrix(rows, cols, upper_bound, lower_bound = 0):
    return [[np.random.uniform(low=lower_bound, high=upper_bound) for _ in range(cols)] for _ in range(rows)]

"""
symnmf_handle - handles the code for goal = "symnmf"
@dataset: the dataset from the cmd input file (already formatted as lists)
@k: (int) number of clusters
@n: (int) length of dataset 
"""
def symnmf_handle(dataset, k, n):
    w_mat = symnmf.norm(dataset)
    w_mat_avg_val = get_average_val_of_matrix(w_mat)
    upper_val = 2 * ((w_mat_avg_val / k)**0.5)
    h_mat = create_random_matrix(n, k, upper_val)

    # call the symnmf() function in module
    res_mat = symnmf.symnmf(h_mat, w_mat, k)

    return res_mat

def get_goal_matrix(goal, k, path):
    dataset = txt_input_to_list(path)
    n = len(dataset)

    if k >= n:
        general_error()
    
    if goal == "symnmf":
        mat = symnmf_handle(dataset, k, n)
    elif goal == "sym":
        mat = symnmf.sym(dataset)
    elif goal == "ddg":
        mat = symnmf.ddg(dataset)
    elif goal == "norm":
        mat = symnmf.norm(dataset)

    return mat
    
def main():
    k, goal, path = set_data(sys.argv)
    dataset = txt_input_to_list(path)
    """
    k, goal = 2, "symnmf"
    dataset = txt_input_to_list("input_1_short.txt")
    """
    # validate goal
    if goal not in ['sym', 'ddg', 'symnmf', 'norm']:
        general_error()

    # validate k
    if (is_float(k)) and (float(k) % 1 != 0):
        general_error()
    
    mat = get_goal_matrix(goal, int(k), path)

    mat_str = matrix_to_str(mat)
    print(mat_str)

if __name__ == "__main__":
    main()
