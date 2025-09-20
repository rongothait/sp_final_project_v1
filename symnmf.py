import numpy as np
import sys
import symnmf_mod as symnmf
import numpy as np

ERR_MSG = "An Error Has Occured"

# set the seed as requested
np.random.seed(1234)

def general_error(msg = ERR_MSG):
    """
    General Error Handeling
    @msg (string, optioanl): the error message. used only for debugging
    """
    print(ERR_MSG)
    exit(1)


def is_float(x):
    """
    is_float - Checks if a string x represents a floating number
    @x: the string
    """
    try:
        float(x)
        return True
    except:
        return False


def get_average_val_of_matrix(mat):
    """
    get_average_val_of_matrix - Calculates the average value of all of the cells in matrix
    @mat: the matrix
    """
    sum = 0
    for row in mat:
        for val in row:
            sum += val
    cnt = len(mat) * len(mat[0])
    return sum / cnt

def validate_str_is_integer(str):
    try:
        k = float(str)
    except:
        general_error("k is not a float : {}".format(str))
    if k % 1 != 0 or k <= 1:
        general_error("k is not an integer number : {}".format(str))
    k = int(k)
    return k


def set_data(data):
    """
    set_data - checks cmd args valididty
    @data: an array of cmd args
    """
    # 1. making sure the input is in correct length
    if len(data) != 4:
        general_error("length of data is incorrect!, expected length of 4 received {}".format(len(data)))

    # validate k
    k = validate_str_is_integer(data[1])
    
    # 2. read arguments from CMD - we can assume they are valid
    goal = data[2]  # 'symnf' / 'sym' / 'ddg' / 'norm'
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        general_error("goal not valid")

    path = data[3]  # path to data set

    return k, goal, path


def txt_input_to_list(path):
    """
    txt_input_to_list - converts the path given to a list of lists
    @path: the path of the file
    """
    try:
        with open(path, 'r') as f:
            lines = [line.strip() for line in f]
            lines = [[float(x) for x in line.split(',')] for line in lines]
    except:
        general_error("problem with opening the file")
    
    return lines


def matrix_to_str(mat):
    """
    matrix_to_str - converts a list of lists (float) to formatted str ready for printing
    @mat: the matrix to convert
    """
    rows = []
    for row in mat:
        row_str = ','.join([f'{x:.4f}' for x in row])
        rows.append(row_str)

    return '\n'.join(rows)


def create_random_matrix(rows, cols, upper_bound, lower_bound = 0):
    """
    Creates a matrix (list of lists) of shape (rows, cols) filled with random float values
    sampled uniformly from the interval [lower_bound, upper_bound).
    @rows (int): number of rows in matrix
    @cols (int): number of cols in matrix
    @upper_bound (float): Upper bound for the random values
    @lower_bound (float): lower bound for the random values
    """
    return [[np.random.uniform(low=lower_bound, high=upper_bound) for _ in range(cols)] for _ in range(rows)]


def init_H(w_mat, k, np_arr = True):
    """
    Initializes the H matrix for SymnNMF using random values based on the average value of w_mat
    w_mat (list): W matrix
    k (int): number of clusters
    np_arr (bool, optional) : if True, returns a NumPy array, otherwise a Python list of lists.
    """
    w_mat_avg_val = get_average_val_of_matrix(w_mat)
    upper_val = 2 * ((w_mat_avg_val / k)**0.5)
    n = len(w_mat)
    init_h_mat =  create_random_matrix(n, k, upper_val)

    if np_arr:
        return np.array(init_h_mat)
    else:
        return init_h_mat


def symnmf_handle(dataset, k, n):
    """
    symnmf_handle - handles the code for goal = "symnmf"
    @dataset: the dataset from the cmd input file (already formatted as lists)
    @k: (int) number of clusters
    @n: (int) length of dataset 
    """
    try:
        w_mat_res = symnmf.norm(dataset)
    except:
        general_error()
    h_mat = init_H(w_mat_res, k, False)

    # call the symnmf() function in module
    try:
        res_mat = symnmf.symnmf(h_mat, w_mat_res, k)
    except:
        general_error()
    return res_mat


def get_goal_matrix(goal, k, dataset):
    """
    handles the call for getting the requested matrix
    @goal (string): the requested goal
    @k (int): number of clusters requested
    @dataset (list): python list of list 
    """
    n = len(dataset)

    if k >= n:  # not allowed
        general_error("k is larget than n")
    
    if goal == "symnmf":
        mat = symnmf_handle(dataset, k, n)
    elif goal == "sym":
        try:
            mat = symnmf.sym(dataset)
        except:
            general_error()
    elif goal == "ddg":
        try:
            mat = symnmf.ddg(dataset)
        except:
            general_error()
    elif goal == "norm":
        try:
            mat = symnmf.norm(dataset)
        except:
            general_error()

    return mat
    
def main():
    k, goal, path = set_data(sys.argv)
    dataset = txt_input_to_list(path)
    
    res = get_goal_matrix(goal, int(k), dataset)

    mat_str = matrix_to_str(res)
    print(mat_str)

if __name__ == "__main__":
    main()
